!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_0d(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     unlim_dim_level, &
                                     corner &

                                    )

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*),  intent(in) :: vdata !< Data that will
                                                      !! be written out
                                                      !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(0), intent(in), optional :: corner !< Array of starting
                                                               !! indices describing
                                                               !! where the data
                                                               !! will be written to.



  integer :: xdim_index
  integer :: ydim_index


  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner &

                         )
    return
  else
    call error("this branch should never be reached.")
  endif

end subroutine domain_write_0d
!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_1d(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     unlim_dim_level, &
                                     corner &
                                    ,edge_lengths &
                                    )

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


  integer :: xdim_index
  integer :: ydim_index


  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner &
                         ,edge_lengths=edge_lengths &
                         )
    return
  else
    call error("this branch should never be reached.")
  endif

end subroutine domain_write_1d
!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_2d(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     unlim_dim_level, &
                                     corner &
                                    ,edge_lengths &
                                    )

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
  integer :: ndims
  integer, dimension(:), allocatable :: dimension_sizes
  integer :: dpos
  integer :: i
  integer :: ds
  integer, dimension(:), allocatable :: isc
  integer, dimension(:), allocatable :: icsize
  integer, dimension(:), allocatable :: jsc
  integer, dimension(:), allocatable :: jcsize
  integer, dimension(2) :: c
  integer, dimension(2) :: e
  integer(kind=int32), dimension(:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:), allocatable :: buf_real64
  integer(kind=int32), dimension(:,:), allocatable :: lbuf_int32
  integer(kind=int64), dimension(:,:), allocatable :: lbuf_int64
  real(kind=real32), dimension(:,:), allocatable :: lbuf_real32
  real(kind=real64), dimension(:,:), allocatable :: lbuf_real64

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner &
                         ,edge_lengths=edge_lengths &
                         )
    return

  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.true.)
  allocate(dimension_sizes(ndims))
  call get_variable_size(fileobj, variable_name, dimension_sizes, broadcast=.true.)
  dpos = get_domain_decomposed_variable_position(fileobj, variable_name)
  do i = 1, size(dimension_sizes)
    if (i .eq. xdim_index) then
      call mpp_get_compute_domain(io_domain, xsize=ds, position=dpos)
      if (size(vdata,i) .ne. ds) then
        call error("domain decomposed x dimension size mismatch.")
      endif
    elseif (i .eq. ydim_index) then
      call mpp_get_compute_domain(io_domain, ysize=ds, position=dpos)
      if (size(vdata,i) .ne. ds) then
        call error("domain decomposed y dimension size mismatch.")
      endif
    endif
  enddo

  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    allocate(isc(size(fileobj%pelist)))
    allocate(icsize(size(fileobj%pelist)))
    allocate(jsc(size(fileobj%pelist)))
    allocate(jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=isc, xsize=icsize, ybegin=jsc, &
                                 ysize=jcsize, position=dpos)
    if (isc(1) .ne. 1) then
      isc(:) = isc(:) - isc(1) + 1
    endif
    if (jsc(1) .ne. 1) then
      jsc(:) = jsc(:) - jsc(1) + 1
    endif
    c(:) = 1
    c(xdim_index) = isc(1)
    c(ydim_index) = jsc(1)
    e(:) = shape(vdata)
    e(xdim_index) = icsize(1)
    e(ydim_index) = jcsize(1)
    select type(vdata)
      type is (integer(kind=int32))
        call allocate_array(buf_int32, dimension_sizes)
        call put_array_section(vdata, buf_int32, c, e)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, dimension_sizes)
        call put_array_section(vdata, buf_int64, c, e)
      type is (real(kind=real32))
        call allocate_array(buf_real32, dimension_sizes)
        call put_array_section(vdata, buf_real32, c, e)
      type is (real(kind=real64))
        call allocate_array(buf_real64, dimension_sizes)
        call put_array_section(vdata, buf_real64, c, e)
      class default
        call error("unsupported type.")
    end select
    do i = 2, size(fileobj%pelist)
      c(xdim_index) = isc(i)
      c(ydim_index) = jsc(i)
      e(xdim_index) = icsize(i)
      e(ydim_index) = jcsize(i)
      select type(vdata)
        type is (integer(kind=int32))
          call allocate_array(lbuf_int32, e)
          call mpp_recv(lbuf_int32, size(lbuf_int32), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_int32, buf_int32, c, e)
          deallocate(lbuf_int32)
        type is (integer(kind=int64))
          call allocate_array(lbuf_int64, e)
          call mpp_recv(lbuf_int64, size(lbuf_int64), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_int64, buf_int64, c, e)
          deallocate(lbuf_int64)
        type is (real(kind=real32))
          call allocate_array(lbuf_real32, e)
          call mpp_recv(lbuf_real32, size(lbuf_real32), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_real32, buf_real32, c, e)
          deallocate(lbuf_real32)
        type is (real(kind=real64))
          call allocate_array(lbuf_real64, e)
          call mpp_recv(lbuf_real64, size(lbuf_real64), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_real64, buf_real64, c, e)
          deallocate(lbuf_real64)
        class default
          call error("unsupported type.")
      end select
    enddo
    select type(vdata)
      type is (integer(kind=int32))
        call netcdf_write_data(fileobj, variable_name, buf_int32, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_int32)
      type is (integer(kind=int64))
        call netcdf_write_data(fileobj, variable_name, buf_int64, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_int64)
      type is (real(kind=real32))
        call netcdf_write_data(fileobj, variable_name, buf_real32, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_real32)
      type is (real(kind=real64))
        call netcdf_write_data(fileobj, variable_name, buf_real64, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_real64)
      class default
        call error("unsupported type.")
    end select
    deallocate(isc)
    deallocate(icsize)
    deallocate(jsc)
    deallocate(jcsize)
  else
    select type(vdata)
      type is (integer(kind=int32))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      type is (integer(kind=int64))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      type is (real(kind=real32))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      type is (real(kind=real64))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      class default
        call error("unsupported type.")
    end select
    call mpp_sync_self(check=EVENT_SEND)
  endif
  deallocate(dimension_sizes)
end subroutine domain_write_2d
!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_3d(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     unlim_dim_level, &
                                     corner &
                                    ,edge_lengths &
                                    )

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
  integer :: ndims
  integer, dimension(:), allocatable :: dimension_sizes
  integer :: dpos
  integer :: i
  integer :: ds
  integer, dimension(:), allocatable :: isc
  integer, dimension(:), allocatable :: icsize
  integer, dimension(:), allocatable :: jsc
  integer, dimension(:), allocatable :: jcsize
  integer, dimension(3) :: c
  integer, dimension(3) :: e
  integer(kind=int32), dimension(:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:), allocatable :: buf_real64
  integer(kind=int32), dimension(:,:,:), allocatable :: lbuf_int32
  integer(kind=int64), dimension(:,:,:), allocatable :: lbuf_int64
  real(kind=real32), dimension(:,:,:), allocatable :: lbuf_real32
  real(kind=real64), dimension(:,:,:), allocatable :: lbuf_real64

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner &
                         ,edge_lengths=edge_lengths &
                         )
    return

  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.true.)
  allocate(dimension_sizes(ndims))
  call get_variable_size(fileobj, variable_name, dimension_sizes, broadcast=.true.)
  dpos = get_domain_decomposed_variable_position(fileobj, variable_name)
  do i = 1, size(dimension_sizes)
    if (i .eq. xdim_index) then
      call mpp_get_compute_domain(io_domain, xsize=ds, position=dpos)
      if (size(vdata,i) .ne. ds) then
        call error("domain decomposed x dimension size mismatch.")
      endif
    elseif (i .eq. ydim_index) then
      call mpp_get_compute_domain(io_domain, ysize=ds, position=dpos)
      if (size(vdata,i) .ne. ds) then
        call error("domain decomposed y dimension size mismatch.")
      endif
    endif
  enddo

  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    allocate(isc(size(fileobj%pelist)))
    allocate(icsize(size(fileobj%pelist)))
    allocate(jsc(size(fileobj%pelist)))
    allocate(jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=isc, xsize=icsize, ybegin=jsc, &
                                 ysize=jcsize, position=dpos)
    if (isc(1) .ne. 1) then
      isc(:) = isc(:) - isc(1) + 1
    endif
    if (jsc(1) .ne. 1) then
      jsc(:) = jsc(:) - jsc(1) + 1
    endif
    c(:) = 1
    c(xdim_index) = isc(1)
    c(ydim_index) = jsc(1)
    e(:) = shape(vdata)
    e(xdim_index) = icsize(1)
    e(ydim_index) = jcsize(1)
    select type(vdata)
      type is (integer(kind=int32))
        call allocate_array(buf_int32, dimension_sizes)
        call put_array_section(vdata, buf_int32, c, e)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, dimension_sizes)
        call put_array_section(vdata, buf_int64, c, e)
      type is (real(kind=real32))
        call allocate_array(buf_real32, dimension_sizes)
        call put_array_section(vdata, buf_real32, c, e)
      type is (real(kind=real64))
        call allocate_array(buf_real64, dimension_sizes)
        call put_array_section(vdata, buf_real64, c, e)
      class default
        call error("unsupported type.")
    end select
    do i = 2, size(fileobj%pelist)
      c(xdim_index) = isc(i)
      c(ydim_index) = jsc(i)
      e(xdim_index) = icsize(i)
      e(ydim_index) = jcsize(i)
      select type(vdata)
        type is (integer(kind=int32))
          call allocate_array(lbuf_int32, e)
          call mpp_recv(lbuf_int32, size(lbuf_int32), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_int32, buf_int32, c, e)
          deallocate(lbuf_int32)
        type is (integer(kind=int64))
          call allocate_array(lbuf_int64, e)
          call mpp_recv(lbuf_int64, size(lbuf_int64), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_int64, buf_int64, c, e)
          deallocate(lbuf_int64)
        type is (real(kind=real32))
          call allocate_array(lbuf_real32, e)
          call mpp_recv(lbuf_real32, size(lbuf_real32), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_real32, buf_real32, c, e)
          deallocate(lbuf_real32)
        type is (real(kind=real64))
          call allocate_array(lbuf_real64, e)
          call mpp_recv(lbuf_real64, size(lbuf_real64), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_real64, buf_real64, c, e)
          deallocate(lbuf_real64)
        class default
          call error("unsupported type.")
      end select
    enddo
    select type(vdata)
      type is (integer(kind=int32))
        call netcdf_write_data(fileobj, variable_name, buf_int32, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_int32)
      type is (integer(kind=int64))
        call netcdf_write_data(fileobj, variable_name, buf_int64, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_int64)
      type is (real(kind=real32))
        call netcdf_write_data(fileobj, variable_name, buf_real32, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_real32)
      type is (real(kind=real64))
        call netcdf_write_data(fileobj, variable_name, buf_real64, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_real64)
      class default
        call error("unsupported type.")
    end select
    deallocate(isc)
    deallocate(icsize)
    deallocate(jsc)
    deallocate(jcsize)
  else
    select type(vdata)
      type is (integer(kind=int32))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      type is (integer(kind=int64))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      type is (real(kind=real32))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      type is (real(kind=real64))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      class default
        call error("unsupported type.")
    end select
    call mpp_sync_self(check=EVENT_SEND)
  endif
  deallocate(dimension_sizes)
end subroutine domain_write_3d
!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_4d(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     unlim_dim_level, &
                                     corner &
                                    ,edge_lengths &
                                    )

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
  integer :: ndims
  integer, dimension(:), allocatable :: dimension_sizes
  integer :: dpos
  integer :: i
  integer :: ds
  integer, dimension(:), allocatable :: isc
  integer, dimension(:), allocatable :: icsize
  integer, dimension(:), allocatable :: jsc
  integer, dimension(:), allocatable :: jcsize
  integer, dimension(4) :: c
  integer, dimension(4) :: e
  integer(kind=int32), dimension(:,:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:,:), allocatable :: buf_real64
  integer(kind=int32), dimension(:,:,:,:), allocatable :: lbuf_int32
  integer(kind=int64), dimension(:,:,:,:), allocatable :: lbuf_int64
  real(kind=real32), dimension(:,:,:,:), allocatable :: lbuf_real32
  real(kind=real64), dimension(:,:,:,:), allocatable :: lbuf_real64

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner &
                         ,edge_lengths=edge_lengths &
                         )
    return

  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.true.)
  allocate(dimension_sizes(ndims))
  call get_variable_size(fileobj, variable_name, dimension_sizes, broadcast=.true.)
  dpos = get_domain_decomposed_variable_position(fileobj, variable_name)
  do i = 1, size(dimension_sizes)
    if (i .eq. xdim_index) then
      call mpp_get_compute_domain(io_domain, xsize=ds, position=dpos)
      if (size(vdata,i) .ne. ds) then
        call error("domain decomposed x dimension size mismatch.")
      endif
    elseif (i .eq. ydim_index) then
      call mpp_get_compute_domain(io_domain, ysize=ds, position=dpos)
      if (size(vdata,i) .ne. ds) then
        call error("domain decomposed y dimension size mismatch.")
      endif
    endif
  enddo

  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    allocate(isc(size(fileobj%pelist)))
    allocate(icsize(size(fileobj%pelist)))
    allocate(jsc(size(fileobj%pelist)))
    allocate(jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=isc, xsize=icsize, ybegin=jsc, &
                                 ysize=jcsize, position=dpos)
    if (isc(1) .ne. 1) then
      isc(:) = isc(:) - isc(1) + 1
    endif
    if (jsc(1) .ne. 1) then
      jsc(:) = jsc(:) - jsc(1) + 1
    endif
    c(:) = 1
    c(xdim_index) = isc(1)
    c(ydim_index) = jsc(1)
    e(:) = shape(vdata)
    e(xdim_index) = icsize(1)
    e(ydim_index) = jcsize(1)
    select type(vdata)
      type is (integer(kind=int32))
        call allocate_array(buf_int32, dimension_sizes)
        call put_array_section(vdata, buf_int32, c, e)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, dimension_sizes)
        call put_array_section(vdata, buf_int64, c, e)
      type is (real(kind=real32))
        call allocate_array(buf_real32, dimension_sizes)
        call put_array_section(vdata, buf_real32, c, e)
      type is (real(kind=real64))
        call allocate_array(buf_real64, dimension_sizes)
        call put_array_section(vdata, buf_real64, c, e)
      class default
        call error("unsupported type.")
    end select
    do i = 2, size(fileobj%pelist)
      c(xdim_index) = isc(i)
      c(ydim_index) = jsc(i)
      e(xdim_index) = icsize(i)
      e(ydim_index) = jcsize(i)
      select type(vdata)
        type is (integer(kind=int32))
          call allocate_array(lbuf_int32, e)
          call mpp_recv(lbuf_int32, size(lbuf_int32), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_int32, buf_int32, c, e)
          deallocate(lbuf_int32)
        type is (integer(kind=int64))
          call allocate_array(lbuf_int64, e)
          call mpp_recv(lbuf_int64, size(lbuf_int64), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_int64, buf_int64, c, e)
          deallocate(lbuf_int64)
        type is (real(kind=real32))
          call allocate_array(lbuf_real32, e)
          call mpp_recv(lbuf_real32, size(lbuf_real32), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_real32, buf_real32, c, e)
          deallocate(lbuf_real32)
        type is (real(kind=real64))
          call allocate_array(lbuf_real64, e)
          call mpp_recv(lbuf_real64, size(lbuf_real64), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_real64, buf_real64, c, e)
          deallocate(lbuf_real64)
        class default
          call error("unsupported type.")
      end select
    enddo
    select type(vdata)
      type is (integer(kind=int32))
        call netcdf_write_data(fileobj, variable_name, buf_int32, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_int32)
      type is (integer(kind=int64))
        call netcdf_write_data(fileobj, variable_name, buf_int64, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_int64)
      type is (real(kind=real32))
        call netcdf_write_data(fileobj, variable_name, buf_real32, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_real32)
      type is (real(kind=real64))
        call netcdf_write_data(fileobj, variable_name, buf_real64, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_real64)
      class default
        call error("unsupported type.")
    end select
    deallocate(isc)
    deallocate(icsize)
    deallocate(jsc)
    deallocate(jcsize)
  else
    select type(vdata)
      type is (integer(kind=int32))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      type is (integer(kind=int64))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      type is (real(kind=real32))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      type is (real(kind=real64))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      class default
        call error("unsupported type.")
    end select
    call mpp_sync_self(check=EVENT_SEND)
  endif
  deallocate(dimension_sizes)
end subroutine domain_write_4d
!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_5d(fileobj, &
                                     variable_name, &
                                     vdata, &
                                     unlim_dim_level, &
                                     corner &
                                    ,edge_lengths &
                                    )

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
  integer :: ndims
  integer, dimension(:), allocatable :: dimension_sizes
  integer :: dpos
  integer :: i
  integer :: ds
  integer, dimension(:), allocatable :: isc
  integer, dimension(:), allocatable :: icsize
  integer, dimension(:), allocatable :: jsc
  integer, dimension(:), allocatable :: jcsize
  integer, dimension(5) :: c
  integer, dimension(5) :: e
  integer(kind=int32), dimension(:,:,:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:,:,:), allocatable :: buf_real64
  integer(kind=int32), dimension(:,:,:,:,:), allocatable :: lbuf_int32
  integer(kind=int64), dimension(:,:,:,:,:), allocatable :: lbuf_int64
  real(kind=real32), dimension(:,:,:,:,:), allocatable :: lbuf_real32
  real(kind=real64), dimension(:,:,:,:,:), allocatable :: lbuf_real64

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner &
                         ,edge_lengths=edge_lengths &
                         )
    return

  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.true.)
  allocate(dimension_sizes(ndims))
  call get_variable_size(fileobj, variable_name, dimension_sizes, broadcast=.true.)
  dpos = get_domain_decomposed_variable_position(fileobj, variable_name)
  do i = 1, size(dimension_sizes)
    if (i .eq. xdim_index) then
      call mpp_get_compute_domain(io_domain, xsize=ds, position=dpos)
      if (size(vdata,i) .ne. ds) then
        call error("domain decomposed x dimension size mismatch.")
      endif
    elseif (i .eq. ydim_index) then
      call mpp_get_compute_domain(io_domain, ysize=ds, position=dpos)
      if (size(vdata,i) .ne. ds) then
        call error("domain decomposed y dimension size mismatch.")
      endif
    endif
  enddo

  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    allocate(isc(size(fileobj%pelist)))
    allocate(icsize(size(fileobj%pelist)))
    allocate(jsc(size(fileobj%pelist)))
    allocate(jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=isc, xsize=icsize, ybegin=jsc, &
                                 ysize=jcsize, position=dpos)
    if (isc(1) .ne. 1) then
      isc(:) = isc(:) - isc(1) + 1
    endif
    if (jsc(1) .ne. 1) then
      jsc(:) = jsc(:) - jsc(1) + 1
    endif
    c(:) = 1
    c(xdim_index) = isc(1)
    c(ydim_index) = jsc(1)
    e(:) = shape(vdata)
    e(xdim_index) = icsize(1)
    e(ydim_index) = jcsize(1)
    select type(vdata)
      type is (integer(kind=int32))
        call allocate_array(buf_int32, dimension_sizes)
        call put_array_section(vdata, buf_int32, c, e)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, dimension_sizes)
        call put_array_section(vdata, buf_int64, c, e)
      type is (real(kind=real32))
        call allocate_array(buf_real32, dimension_sizes)
        call put_array_section(vdata, buf_real32, c, e)
      type is (real(kind=real64))
        call allocate_array(buf_real64, dimension_sizes)
        call put_array_section(vdata, buf_real64, c, e)
      class default
        call error("unsupported type.")
    end select
    do i = 2, size(fileobj%pelist)
      c(xdim_index) = isc(i)
      c(ydim_index) = jsc(i)
      e(xdim_index) = icsize(i)
      e(ydim_index) = jcsize(i)
      select type(vdata)
        type is (integer(kind=int32))
          call allocate_array(lbuf_int32, e)
          call mpp_recv(lbuf_int32, size(lbuf_int32), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_int32, buf_int32, c, e)
          deallocate(lbuf_int32)
        type is (integer(kind=int64))
          call allocate_array(lbuf_int64, e)
          call mpp_recv(lbuf_int64, size(lbuf_int64), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_int64, buf_int64, c, e)
          deallocate(lbuf_int64)
        type is (real(kind=real32))
          call allocate_array(lbuf_real32, e)
          call mpp_recv(lbuf_real32, size(lbuf_real32), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_real32, buf_real32, c, e)
          deallocate(lbuf_real32)
        type is (real(kind=real64))
          call allocate_array(lbuf_real64, e)
          call mpp_recv(lbuf_real64, size(lbuf_real64), fileobj%pelist(i))
          call mpp_sync_self(check=EVENT_RECV)
          call put_array_section(lbuf_real64, buf_real64, c, e)
          deallocate(lbuf_real64)
        class default
          call error("unsupported type.")
      end select
    enddo
    select type(vdata)
      type is (integer(kind=int32))
        call netcdf_write_data(fileobj, variable_name, buf_int32, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_int32)
      type is (integer(kind=int64))
        call netcdf_write_data(fileobj, variable_name, buf_int64, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_int64)
      type is (real(kind=real32))
        call netcdf_write_data(fileobj, variable_name, buf_real32, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_real32)
      type is (real(kind=real64))
        call netcdf_write_data(fileobj, variable_name, buf_real64, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(buf_real64)
      class default
        call error("unsupported type.")
    end select
    deallocate(isc)
    deallocate(icsize)
    deallocate(jsc)
    deallocate(jcsize)
  else
    select type(vdata)
      type is (integer(kind=int32))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      type is (integer(kind=int64))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      type is (real(kind=real32))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      type is (real(kind=real64))
        call mpp_send(vdata, size(vdata), fileobj%io_root)
      class default
        call error("unsupported type.")
    end select
    call mpp_sync_self(check=EVENT_SEND)
  endif
  deallocate(dimension_sizes)
end subroutine domain_write_5d
