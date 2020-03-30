!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_0d(fileobj, variable_name, cdata, unlim_dim_level, &
                               corner)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), intent(in) :: cdata !< Compressed data that
                                !! will be gathered and
                                !! written to the
                                !! netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, intent(in), optional :: corner !< Array of starting
                                          !! indices describing
                                          !! where the data
                                          !! will be written to.

  integer, dimension(2) :: compressed_dim_index

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)
  if (compressed_dim_index(1) .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner)
    return
  else
    call error("this branch should never be reached.")
  endif
end subroutine compressed_write_0d


!> @brief Wrapper to distinguish interfaces.
subroutine compressed_write_0d_wrap(fileobj, variable_name, cdata, unlim_dim_level, &
                                    corner)

  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), intent(in) :: cdata !< Compressed data that
                                !! will be gathered and
                                !! written to the
                                !! netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, intent(in), optional :: corner !< Array of starting
                                          !! indices describing
                                          !! where the data
                                          !! will be written to.

  call compressed_write(fileobj, variable_name, cdata, unlim_dim_level, corner)
end subroutine compressed_write_0d_wrap


!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_1d(fileobj, variable_name, cdata, unlim_dim_level, &
                               corner, edge_lengths)

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

  integer, dimension(2) :: compressed_dim_index
  integer, dimension(1) :: c
  integer, dimension(1) :: e
  integer :: i
  integer(kind=int32), dimension(:), allocatable :: buf_int32
  integer(kind=int64), dimension(:), allocatable :: buf_int64
  real(kind=real32), dimension(:), allocatable :: buf_real32
  real(kind=real64), dimension(:), allocatable :: buf_real64

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)

  if (compressed_dim_index(1) .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner, &
                           edge_lengths=edge_lengths)
    return
  endif

  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    c(:) = 1
    e(:) = shape(cdata)
    do i = 1, size(fileobj%pelist)
      c(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_corner(i)
      e(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems(i)
      if (i .eq. 1) then
        call netcdf_write_data(fileobj, variable_name, cdata, &
                               unlim_dim_level=unlim_dim_level, corner=c, &
                               edge_lengths=e)
      else
        select type(cdata)
          type is (integer(kind=int32))
            call allocate_array(buf_int32, e)
            call mpp_recv(buf_int32, size(buf_int32), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_int32, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_int32)
          type is (integer(kind=int64))
            call allocate_array(buf_int64, e)
            call mpp_recv(buf_int64, size(buf_int64), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_int64, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_int64)
          type is (real(kind=real32))
            call allocate_array(buf_real32, e)
            call mpp_recv(buf_real32, size(buf_real32), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_real32, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_real32)
          type is (real(kind=real64))
            call allocate_array(buf_real64, e)
            call mpp_recv(buf_real64, size(buf_real64), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_real64, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_real64)
          class default
            call error("unsupported type.")
        end select
      endif
    enddo
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


!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_2d(fileobj, variable_name, cdata, unlim_dim_level, &
                               corner, edge_lengths)

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

  integer, dimension(2) :: compressed_dim_index
  integer, dimension(2) :: c
  integer, dimension(2) :: e
  integer :: i
  integer(kind=int32), dimension(:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:), allocatable :: buf_real64

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)
  if (compressed_dim_index(1) .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner, &
                           edge_lengths=edge_lengths)
    return
  endif

  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    c(:) = 1
    e(:) = shape(cdata)
    do i = 1, size(fileobj%pelist)
      c(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_corner(i)
      e(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems(i)
      if (i .eq. 1) then
        call netcdf_write_data(fileobj, variable_name, cdata, &
                               unlim_dim_level=unlim_dim_level, corner=c, &
                               edge_lengths=e)
      else
        select type(cdata)
          type is (integer(kind=int32))
            call allocate_array(buf_int32, e)
            call mpp_recv(buf_int32, size(buf_int32), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_int32, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_int32)
          type is (integer(kind=int64))
            call allocate_array(buf_int64, e)
            call mpp_recv(buf_int64, size(buf_int64), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_int64, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_int64)
          type is (real(kind=real32))
            call allocate_array(buf_real32, e)
            call mpp_recv(buf_real32, size(buf_real32), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_real32, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_real32)
          type is (real(kind=real64))
            call allocate_array(buf_real64, e)
            call mpp_recv(buf_real64, size(buf_real64), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_real64, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_real64)
          class default
            call error("unsupported type.")
        end select
      endif
    enddo
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


!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_3d(fileobj, variable_name, cdata, unlim_dim_level, &
                               corner, edge_lengths)

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

  integer, dimension(2) :: compressed_dim_index
  integer, dimension(3) :: c
  integer, dimension(3) :: e
  integer :: i
  integer(kind=int32), dimension(:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:), allocatable :: buf_real64

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)
  if (compressed_dim_index(1) .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner, &
                           edge_lengths=edge_lengths)
    return
  endif

  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    c(:) = 1
    e(:) = shape(cdata)
    do i = 1, size(fileobj%pelist)
      c(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_corner(i)
      e(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems(i)
      if (i .eq. 1) then
        call netcdf_write_data(fileobj, variable_name, cdata, &
                               unlim_dim_level=unlim_dim_level, corner=c, &
                               edge_lengths=e)
      else
        select type(cdata)
          type is (integer(kind=int32))
            call allocate_array(buf_int32, e)
            call mpp_recv(buf_int32, size(buf_int32), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_int32, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_int32)
          type is (integer(kind=int64))
            call allocate_array(buf_int64, e)
            call mpp_recv(buf_int64, size(buf_int64), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_int64, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_int64)
          type is (real(kind=real32))
            call allocate_array(buf_real32, e)
            call mpp_recv(buf_real32, size(buf_real32), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_real32, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_real32)
          type is (real(kind=real64))
            call allocate_array(buf_real64, e)
            call mpp_recv(buf_real64, size(buf_real64), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_real64, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_real64)
          class default
            call error("unsupported type.")
        end select
      endif
    enddo
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


!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_4d(fileobj, variable_name, cdata, unlim_dim_level, &
                               corner, edge_lengths)

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

  integer, dimension(2) :: compressed_dim_index
  integer, dimension(4) :: c
  integer, dimension(4) :: e
  integer :: i
  integer(kind=int32), dimension(:,:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:,:), allocatable :: buf_real64

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)
  if (compressed_dim_index(1) .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner, &
                           edge_lengths=edge_lengths)
    return
  endif

  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    c(:) = 1
    e(:) = shape(cdata)
    do i = 1, size(fileobj%pelist)
      c(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_corner(i)
      e(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems(i)
      if (i .eq. 1) then
        call netcdf_write_data(fileobj, variable_name, cdata, &
                               unlim_dim_level=unlim_dim_level, corner=c, &
                               edge_lengths=e)
      else
        select type(cdata)
          type is (integer(kind=int32))
            call allocate_array(buf_int32, e)
            call mpp_recv(buf_int32, size(buf_int32), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_int32, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_int32)
          type is (integer(kind=int64))
            call allocate_array(buf_int64, e)
            call mpp_recv(buf_int64, size(buf_int64), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_int64, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_int64)
          type is (real(kind=real32))
            call allocate_array(buf_real32, e)
            call mpp_recv(buf_real32, size(buf_real32), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_real32, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_real32)
          type is (real(kind=real64))
            call allocate_array(buf_real64, e)
            call mpp_recv(buf_real64, size(buf_real64), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_real64, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_real64)
          class default
            call error("unsupported type.")
        end select
      endif
    enddo
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


!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_5d(fileobj, variable_name, cdata, unlim_dim_level, &
                               corner, edge_lengths)

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

  integer, dimension(2) :: compressed_dim_index
  integer, dimension(5) :: c
  integer, dimension(5) :: e
  integer :: i
  integer(kind=int32), dimension(:,:,:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:,:,:), allocatable :: buf_real64

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)
  if (compressed_dim_index(1) .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner, &
                           edge_lengths=edge_lengths)
    return
  endif

  !Gather the data onto the I/O root and write it out.
  if (fileobj%is_root) then
    c(:) = 1
    e(:) = shape(cdata)
    do i = 1, size(fileobj%pelist)
      c(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_corner(i)
      e(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems(i)
      if (i .eq. 1) then
        call netcdf_write_data(fileobj, variable_name, cdata, &
                               unlim_dim_level=unlim_dim_level, corner=c, &
                               edge_lengths=e)
      else
        select type(cdata)
          type is (integer(kind=int32))
            call allocate_array(buf_int32, e)
            call mpp_recv(buf_int32, size(buf_int32), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_int32, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_int32)
          type is (integer(kind=int64))
            call allocate_array(buf_int64, e)
            call mpp_recv(buf_int64, size(buf_int64), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_int64, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_int64)
          type is (real(kind=real32))
            call allocate_array(buf_real32, e)
            call mpp_recv(buf_real32, size(buf_real32), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_real32, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_real32)
          type is (real(kind=real64))
            call allocate_array(buf_real64, e)
            call mpp_recv(buf_real64, size(buf_real64), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_real64, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_real64)
          class default
            call error("unsupported type.")
        end select
      endif
    enddo
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
subroutine compressed_write_1d_wrap(fileobj, variable_name, cdata, unlim_dim_level, &
                                    corner, edge_lengths)

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

  call compressed_write(fileobj, variable_name, cdata, unlim_dim_level, corner, &
                        edge_lengths=edge_lengths)

end subroutine compressed_write_1d_wrap


!> @brief Wrapper to distinguish interfaces.
subroutine compressed_write_2d_wrap(fileobj, variable_name, cdata, unlim_dim_level, &
                                    corner, edge_lengths)

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

  call compressed_write(fileobj, variable_name, cdata, unlim_dim_level, corner, &
                        edge_lengths=edge_lengths)
end subroutine compressed_write_2d_wrap


!> @brief Wrapper to distinguish interfaces.
subroutine compressed_write_3d_wrap(fileobj, variable_name, cdata, unlim_dim_level, &
                                    corner, edge_lengths)

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

  call compressed_write(fileobj, variable_name, cdata, unlim_dim_level, corner, &
                        edge_lengths=edge_lengths)
end subroutine compressed_write_3d_wrap


!> @brief Wrapper to distinguish interfaces.
subroutine compressed_write_4d_wrap(fileobj, variable_name, cdata, unlim_dim_level, &
                                    corner, edge_lengths)

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

  call compressed_write(fileobj, variable_name, cdata, unlim_dim_level, corner, &
                        edge_lengths=edge_lengths)
end subroutine compressed_write_4d_wrap


!> @brief Wrapper to distinguish interfaces.
subroutine compressed_write_5d_wrap(fileobj, variable_name, cdata, unlim_dim_level, &
                                    corner, edge_lengths)

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

  call compressed_write(fileobj, variable_name, cdata, unlim_dim_level, corner, &
                        edge_lengths=edge_lengths)
end subroutine compressed_write_5d_wrap
