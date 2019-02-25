!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_0d(fileobj, &
                                         variable_name, &
                                         cdata, &
                                         unlim_dim_level, &
                                         corner &

                                        )

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*),  intent(in) :: cdata !< Compressed data that
                                                      !! will be gathered and
                                                      !! written to the 
                                                      !! netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(0), intent(in), optional :: corner !< Array of starting
                                                               !! indices describing
                                                               !! where the data
                                                               !! will be written to.


  integer :: compressed_dim_index


  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)
  if (compressed_dim_index .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner &

                          )
    return
  else
    call error("this branch should never be reached.")
  endif

end subroutine compressed_write_0d


!> @brief Wrapper to distinguish interfaces.
subroutine compressed_write_0d_wrap(fileobj, &
                                              variable_name, &
                                              cdata, &
                                              unlim_dim_level, &
                                              corner &

                                             )

  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*),  intent(in) :: cdata !< Compressed data that
                                                      !! will be gathered and
                                                      !! written to the 
                                                      !! netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(0), intent(in), optional :: corner !< Array of starting
                                                               !! indices describing
                                                               !! where the data
                                                               !! will be written to.


  call compressed_write(fileobj, variable_name, cdata, unlim_dim_level, corner &

                       )
end subroutine compressed_write_0d_wrap
!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_1d(fileobj, &
                                         variable_name, &
                                         cdata, &
                                         unlim_dim_level, &
                                         corner &
                                        ,edge_lengths &
                                        )

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(in) :: cdata !< Compressed data that
                                                      !! will be gathered and
                                                      !! written to the 
                                                      !! netcdf file.
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

  integer :: compressed_dim_index
  integer, dimension(:), allocatable :: dimension_sizes
  integer :: ndims
  character(len=nf90_max_name), dimension(:), allocatable :: dimension_names
  integer :: cindex
  integer :: i
  integer(kind=int32), dimension(:), allocatable :: buf_int32
  integer(kind=int64), dimension(:), allocatable :: buf_int64
  real(kind=real32), dimension(:), allocatable :: buf_real32
  real(kind=real64), dimension(:), allocatable :: buf_real64
  integer(kind=int32), dimension(:), allocatable :: lbuf_int32
  integer(kind=int64), dimension(:), allocatable :: lbuf_int64
  real(kind=real32), dimension(:), allocatable :: lbuf_real32
  real(kind=real64), dimension(:), allocatable :: lbuf_real64
  integer, dimension(1) :: e
  integer, dimension(1) :: c

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)
  if (compressed_dim_index .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner &
                          ,edge_lengths=edge_lengths &
                          )
    return

  endif
  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.false.)
    allocate(dimension_sizes(ndims))
    allocate(dimension_names(ndims))
    call get_variable_size(fileobj, variable_name, dimension_sizes, broadcast=.false.)
    call get_variable_dimension_names(fileobj, variable_name, dimension_names, broadcast=.false.)
    do i = 1, size(dimension_names)
      cindex = get_compressed_dimension_index(fileobj, dimension_names(i))
      if (cindex .ne. dimension_not_found) then
        exit
      endif
    enddo
    c(:) = 1
    c(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_corner(1)
    e(:) = shape(cdata)
    e(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_nelems(1)
    select type(cdata)
      type is (integer(kind=int32))
        call allocate_array(buf_int32, dimension_sizes)
        call put_array_section(cdata, buf_int32, c, e)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, dimension_sizes)
        call put_array_section(cdata, buf_int64, c, e)
      type is (real(kind=real32))
        call allocate_array(buf_real32, dimension_sizes)
        call put_array_section(cdata, buf_real32, c, e)
      type is (real(kind=real64))
        call allocate_array(buf_real64, dimension_sizes)
        call put_array_section(cdata, buf_real64, c, e)
      class default
        call error("unsupported type.")
    end select
    do i = 2, size(fileobj%pelist)
      c(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_corner(i)
      e(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_nelems(i)
      select type(cdata)
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
    select type(cdata)
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
    deallocate(dimension_sizes)
    deallocate(dimension_names)
  else
    select type(cdata)
      type is (integer(kind=int32))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (integer(kind=int64))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=real32))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=real64))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      class default
        call error("unsupported type.")
    end select
    call mpp_sync_self(check=EVENT_SEND)
  endif
end subroutine compressed_write_1d


!> @brief Wrapper to distinguish interfaces.
subroutine compressed_write_1d_wrap(fileobj, &
                                              variable_name, &
                                              cdata, &
                                              unlim_dim_level, &
                                              corner &
                                            ,edge_lengths &
                                             )

  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(in) :: cdata !< Compressed data that
                                                      !! will be gathered and
                                                      !! written to the 
                                                      !! netcdf file.
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

  call compressed_write(fileobj, variable_name, cdata, unlim_dim_level, corner &
                       ,edge_lengths=edge_lengths &
                       )
end subroutine compressed_write_1d_wrap
!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_2d(fileobj, &
                                         variable_name, &
                                         cdata, &
                                         unlim_dim_level, &
                                         corner &
                                        ,edge_lengths &
                                        )

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in) :: cdata !< Compressed data that
                                                      !! will be gathered and
                                                      !! written to the 
                                                      !! netcdf file.
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

  integer :: compressed_dim_index
  integer, dimension(:), allocatable :: dimension_sizes
  integer :: ndims
  character(len=nf90_max_name), dimension(:), allocatable :: dimension_names
  integer :: cindex
  integer :: i
  integer(kind=int32), dimension(:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:), allocatable :: buf_real64
  integer(kind=int32), dimension(:,:), allocatable :: lbuf_int32
  integer(kind=int64), dimension(:,:), allocatable :: lbuf_int64
  real(kind=real32), dimension(:,:), allocatable :: lbuf_real32
  real(kind=real64), dimension(:,:), allocatable :: lbuf_real64
  integer, dimension(2) :: e
  integer, dimension(2) :: c

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)
  if (compressed_dim_index .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner &
                          ,edge_lengths=edge_lengths &
                          )
    return

  endif
  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.false.)
    allocate(dimension_sizes(ndims))
    allocate(dimension_names(ndims))
    call get_variable_size(fileobj, variable_name, dimension_sizes, broadcast=.false.)
    call get_variable_dimension_names(fileobj, variable_name, dimension_names, broadcast=.false.)
    do i = 1, size(dimension_names)
      cindex = get_compressed_dimension_index(fileobj, dimension_names(i))
      if (cindex .ne. dimension_not_found) then
        exit
      endif
    enddo
    c(:) = 1
    c(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_corner(1)
    e(:) = shape(cdata)
    e(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_nelems(1)
    select type(cdata)
      type is (integer(kind=int32))
        call allocate_array(buf_int32, dimension_sizes)
        call put_array_section(cdata, buf_int32, c, e)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, dimension_sizes)
        call put_array_section(cdata, buf_int64, c, e)
      type is (real(kind=real32))
        call allocate_array(buf_real32, dimension_sizes)
        call put_array_section(cdata, buf_real32, c, e)
      type is (real(kind=real64))
        call allocate_array(buf_real64, dimension_sizes)
        call put_array_section(cdata, buf_real64, c, e)
      class default
        call error("unsupported type.")
    end select
    do i = 2, size(fileobj%pelist)
      c(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_corner(i)
      e(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_nelems(i)
      select type(cdata)
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
    select type(cdata)
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
    deallocate(dimension_sizes)
    deallocate(dimension_names)
  else
    select type(cdata)
      type is (integer(kind=int32))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (integer(kind=int64))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=real32))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=real64))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      class default
        call error("unsupported type.")
    end select
    call mpp_sync_self(check=EVENT_SEND)
  endif
end subroutine compressed_write_2d


!> @brief Wrapper to distinguish interfaces.
subroutine compressed_write_2d_wrap(fileobj, &
                                              variable_name, &
                                              cdata, &
                                              unlim_dim_level, &
                                              corner &
                                            ,edge_lengths &
                                             )

  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in) :: cdata !< Compressed data that
                                                      !! will be gathered and
                                                      !! written to the 
                                                      !! netcdf file.
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

  call compressed_write(fileobj, variable_name, cdata, unlim_dim_level, corner &
                       ,edge_lengths=edge_lengths &
                       )
end subroutine compressed_write_2d_wrap
!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_3d(fileobj, &
                                         variable_name, &
                                         cdata, &
                                         unlim_dim_level, &
                                         corner &
                                        ,edge_lengths &
                                        )

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in) :: cdata !< Compressed data that
                                                      !! will be gathered and
                                                      !! written to the 
                                                      !! netcdf file.
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

  integer :: compressed_dim_index
  integer, dimension(:), allocatable :: dimension_sizes
  integer :: ndims
  character(len=nf90_max_name), dimension(:), allocatable :: dimension_names
  integer :: cindex
  integer :: i
  integer(kind=int32), dimension(:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:), allocatable :: buf_real64
  integer(kind=int32), dimension(:,:,:), allocatable :: lbuf_int32
  integer(kind=int64), dimension(:,:,:), allocatable :: lbuf_int64
  real(kind=real32), dimension(:,:,:), allocatable :: lbuf_real32
  real(kind=real64), dimension(:,:,:), allocatable :: lbuf_real64
  integer, dimension(3) :: e
  integer, dimension(3) :: c

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)
  if (compressed_dim_index .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner &
                          ,edge_lengths=edge_lengths &
                          )
    return

  endif
  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.false.)
    allocate(dimension_sizes(ndims))
    allocate(dimension_names(ndims))
    call get_variable_size(fileobj, variable_name, dimension_sizes, broadcast=.false.)
    call get_variable_dimension_names(fileobj, variable_name, dimension_names, broadcast=.false.)
    do i = 1, size(dimension_names)
      cindex = get_compressed_dimension_index(fileobj, dimension_names(i))
      if (cindex .ne. dimension_not_found) then
        exit
      endif
    enddo
    c(:) = 1
    c(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_corner(1)
    e(:) = shape(cdata)
    e(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_nelems(1)
    select type(cdata)
      type is (integer(kind=int32))
        call allocate_array(buf_int32, dimension_sizes)
        call put_array_section(cdata, buf_int32, c, e)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, dimension_sizes)
        call put_array_section(cdata, buf_int64, c, e)
      type is (real(kind=real32))
        call allocate_array(buf_real32, dimension_sizes)
        call put_array_section(cdata, buf_real32, c, e)
      type is (real(kind=real64))
        call allocate_array(buf_real64, dimension_sizes)
        call put_array_section(cdata, buf_real64, c, e)
      class default
        call error("unsupported type.")
    end select
    do i = 2, size(fileobj%pelist)
      c(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_corner(i)
      e(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_nelems(i)
      select type(cdata)
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
    select type(cdata)
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
    deallocate(dimension_sizes)
    deallocate(dimension_names)
  else
    select type(cdata)
      type is (integer(kind=int32))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (integer(kind=int64))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=real32))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=real64))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      class default
        call error("unsupported type.")
    end select
    call mpp_sync_self(check=EVENT_SEND)
  endif
end subroutine compressed_write_3d


!> @brief Wrapper to distinguish interfaces.
subroutine compressed_write_3d_wrap(fileobj, &
                                              variable_name, &
                                              cdata, &
                                              unlim_dim_level, &
                                              corner &
                                            ,edge_lengths &
                                             )

  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in) :: cdata !< Compressed data that
                                                      !! will be gathered and
                                                      !! written to the 
                                                      !! netcdf file.
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

  call compressed_write(fileobj, variable_name, cdata, unlim_dim_level, corner &
                       ,edge_lengths=edge_lengths &
                       )
end subroutine compressed_write_3d_wrap
!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_4d(fileobj, &
                                         variable_name, &
                                         cdata, &
                                         unlim_dim_level, &
                                         corner &
                                        ,edge_lengths &
                                        )

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(in) :: cdata !< Compressed data that
                                                      !! will be gathered and
                                                      !! written to the 
                                                      !! netcdf file.
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

  integer :: compressed_dim_index
  integer, dimension(:), allocatable :: dimension_sizes
  integer :: ndims
  character(len=nf90_max_name), dimension(:), allocatable :: dimension_names
  integer :: cindex
  integer :: i
  integer(kind=int32), dimension(:,:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:,:), allocatable :: buf_real64
  integer(kind=int32), dimension(:,:,:,:), allocatable :: lbuf_int32
  integer(kind=int64), dimension(:,:,:,:), allocatable :: lbuf_int64
  real(kind=real32), dimension(:,:,:,:), allocatable :: lbuf_real32
  real(kind=real64), dimension(:,:,:,:), allocatable :: lbuf_real64
  integer, dimension(4) :: e
  integer, dimension(4) :: c

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)
  if (compressed_dim_index .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner &
                          ,edge_lengths=edge_lengths &
                          )
    return

  endif
  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.false.)
    allocate(dimension_sizes(ndims))
    allocate(dimension_names(ndims))
    call get_variable_size(fileobj, variable_name, dimension_sizes, broadcast=.false.)
    call get_variable_dimension_names(fileobj, variable_name, dimension_names, broadcast=.false.)
    do i = 1, size(dimension_names)
      cindex = get_compressed_dimension_index(fileobj, dimension_names(i))
      if (cindex .ne. dimension_not_found) then
        exit
      endif
    enddo
    c(:) = 1
    c(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_corner(1)
    e(:) = shape(cdata)
    e(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_nelems(1)
    select type(cdata)
      type is (integer(kind=int32))
        call allocate_array(buf_int32, dimension_sizes)
        call put_array_section(cdata, buf_int32, c, e)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, dimension_sizes)
        call put_array_section(cdata, buf_int64, c, e)
      type is (real(kind=real32))
        call allocate_array(buf_real32, dimension_sizes)
        call put_array_section(cdata, buf_real32, c, e)
      type is (real(kind=real64))
        call allocate_array(buf_real64, dimension_sizes)
        call put_array_section(cdata, buf_real64, c, e)
      class default
        call error("unsupported type.")
    end select
    do i = 2, size(fileobj%pelist)
      c(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_corner(i)
      e(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_nelems(i)
      select type(cdata)
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
    select type(cdata)
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
    deallocate(dimension_sizes)
    deallocate(dimension_names)
  else
    select type(cdata)
      type is (integer(kind=int32))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (integer(kind=int64))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=real32))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=real64))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      class default
        call error("unsupported type.")
    end select
    call mpp_sync_self(check=EVENT_SEND)
  endif
end subroutine compressed_write_4d


!> @brief Wrapper to distinguish interfaces.
subroutine compressed_write_4d_wrap(fileobj, &
                                              variable_name, &
                                              cdata, &
                                              unlim_dim_level, &
                                              corner &
                                            ,edge_lengths &
                                             )

  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(in) :: cdata !< Compressed data that
                                                      !! will be gathered and
                                                      !! written to the 
                                                      !! netcdf file.
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

  call compressed_write(fileobj, variable_name, cdata, unlim_dim_level, corner &
                       ,edge_lengths=edge_lengths &
                       )
end subroutine compressed_write_4d_wrap
!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_5d(fileobj, &
                                         variable_name, &
                                         cdata, &
                                         unlim_dim_level, &
                                         corner &
                                        ,edge_lengths &
                                        )

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(in) :: cdata !< Compressed data that
                                                      !! will be gathered and
                                                      !! written to the 
                                                      !! netcdf file.
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

  integer :: compressed_dim_index
  integer, dimension(:), allocatable :: dimension_sizes
  integer :: ndims
  character(len=nf90_max_name), dimension(:), allocatable :: dimension_names
  integer :: cindex
  integer :: i
  integer(kind=int32), dimension(:,:,:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:,:,:), allocatable :: buf_real64
  integer(kind=int32), dimension(:,:,:,:,:), allocatable :: lbuf_int32
  integer(kind=int64), dimension(:,:,:,:,:), allocatable :: lbuf_int64
  real(kind=real32), dimension(:,:,:,:,:), allocatable :: lbuf_real32
  real(kind=real64), dimension(:,:,:,:,:), allocatable :: lbuf_real64
  integer, dimension(5) :: e
  integer, dimension(5) :: c

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)
  if (compressed_dim_index .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner &
                          ,edge_lengths=edge_lengths &
                          )
    return

  endif
  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.false.)
    allocate(dimension_sizes(ndims))
    allocate(dimension_names(ndims))
    call get_variable_size(fileobj, variable_name, dimension_sizes, broadcast=.false.)
    call get_variable_dimension_names(fileobj, variable_name, dimension_names, broadcast=.false.)
    do i = 1, size(dimension_names)
      cindex = get_compressed_dimension_index(fileobj, dimension_names(i))
      if (cindex .ne. dimension_not_found) then
        exit
      endif
    enddo
    c(:) = 1
    c(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_corner(1)
    e(:) = shape(cdata)
    e(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_nelems(1)
    select type(cdata)
      type is (integer(kind=int32))
        call allocate_array(buf_int32, dimension_sizes)
        call put_array_section(cdata, buf_int32, c, e)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, dimension_sizes)
        call put_array_section(cdata, buf_int64, c, e)
      type is (real(kind=real32))
        call allocate_array(buf_real32, dimension_sizes)
        call put_array_section(cdata, buf_real32, c, e)
      type is (real(kind=real64))
        call allocate_array(buf_real64, dimension_sizes)
        call put_array_section(cdata, buf_real64, c, e)
      class default
        call error("unsupported type.")
    end select
    do i = 2, size(fileobj%pelist)
      c(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_corner(i)
      e(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_nelems(i)
      select type(cdata)
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
    select type(cdata)
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
    deallocate(dimension_sizes)
    deallocate(dimension_names)
  else
    select type(cdata)
      type is (integer(kind=int32))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (integer(kind=int64))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=real32))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=real64))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      class default
        call error("unsupported type.")
    end select
    call mpp_sync_self(check=EVENT_SEND)
  endif
end subroutine compressed_write_5d


!> @brief Wrapper to distinguish interfaces.
subroutine compressed_write_5d_wrap(fileobj, &
                                              variable_name, &
                                              cdata, &
                                              unlim_dim_level, &
                                              corner &
                                            ,edge_lengths &
                                             )

  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(in) :: cdata !< Compressed data that
                                                      !! will be gathered and
                                                      !! written to the 
                                                      !! netcdf file.
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

  call compressed_write(fileobj, variable_name, cdata, unlim_dim_level, corner &
                       ,edge_lengths=edge_lengths &
                       )
end subroutine compressed_write_5d_wrap
