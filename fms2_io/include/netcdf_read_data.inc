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

!> @brief Read in data from a variable in a netcdf file.
subroutine netcdf_read_data_0d(fileobj, variable_name, buf, unlim_dim_level, &
                               corner, broadcast)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*),  intent(inout) :: buf !< Array that the data will be read into.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.
  integer, intent(in), optional :: corner !< Index of the string to read if the variable
                                          !! contains a 1D array of strings.
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.

  logical :: bcast
  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer, dimension(1) :: c
  integer, dimension(2) :: start
  integer :: ndims
  character, dimension(:), allocatable :: charbuf
  integer, dimension(:), allocatable :: dimsizes
  integer :: i
  character(len=1024), dimension(1) :: buf1d

  if (present(broadcast)) then
    bcast = broadcast
  else
    bcast = .true.
  endif
  c(1) = 1
  if (present(unlim_dim_level)) then
    unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                             broadcast=bcast)
    if (unlim_dim_index .ne. 1) then
      call error("unlimited dimension must be the slowest varying dimension.")
    endif
    c(unlim_dim_index) = unlim_dim_level
  endif
  if (fileobj%is_root) then
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(buf)
      type is (integer(kind=int32))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c)
      type is (integer(kind=int64))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c)
      type is (real(kind=real32))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c)
      type is (real(kind=real64))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c)
      type is (character(len=*))
        start(:) = 1
        ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.false.)
        allocate(dimsizes(ndims))
        call get_variable_size(fileobj, variable_name, dimsizes, broadcast=.false.)
        allocate(charbuf(dimsizes(1)))
        charbuf(:) = ""
        if (ndims .eq. 2) then
          if (present(corner)) then
            start(2) = corner
          endif
          dimsizes(2) = 1
        elseif (ndims .gt. 2) then
          call error("Only scalar and 1d string values are currently supported.")
        endif
        err = nf90_get_var(fileobj%ncid, varid, charbuf, start=start, count=dimsizes)
        if (len(buf) .lt. dimsizes(1)) then
          call error("character buffer is too small; increase length.")
        endif
        buf = ""
        do i = 1, dimsizes(1)
          if (charbuf(i) .eq. char(0)) then
            exit
          endif
          buf(i:i) = charbuf(i)
        enddo
        deallocate(charbuf)
        deallocate(dimsizes)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)
  endif
  if (bcast) then
    select type(buf)
      type is (integer(kind=int32))
        call mpp_broadcast(buf, fileobj%io_root, pelist=fileobj%pelist)
      type is (integer(kind=int64))
        call mpp_broadcast(buf, fileobj%io_root, pelist=fileobj%pelist)
      type is (real(kind=real32))
        call mpp_broadcast(buf, fileobj%io_root, pelist=fileobj%pelist)
      type is (real(kind=real64))
        call mpp_broadcast(buf, fileobj%io_root, pelist=fileobj%pelist)
      type is (character(len=*))
        call string_copy(buf1d(1), buf)
        call mpp_broadcast(buf1d, len(buf1d(1)), fileobj%io_root, pelist=fileobj%pelist)
        call string_copy(buf, buf1d(1))
      class default
        call error("unsupported type.")
    end select
  endif
end subroutine netcdf_read_data_0d


!> @brief Read in data from a variable in a netcdf file.
subroutine netcdf_read_data_1d(fileobj, variable_name, buf, unlim_dim_level, &
                               corner, edge_lengths, broadcast)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(inout) :: buf !< Array that the data
                                               !! will be read into.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be read from.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be read
                                                              !! in each dimension.
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.

  logical :: bcast
  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer, dimension(2) :: c
  integer, dimension(2) :: e
  integer :: ndims
  integer, dimension(2) :: start
  integer, dimension(2) :: dimsizes
  character, dimension(:,:), allocatable :: charbuf
  character(len=1024) :: sbuf
  integer :: i
  integer :: j

  if (present(broadcast)) then
    bcast = broadcast
  else
    bcast = .true.
  endif
  c(:) = 1
  if (present(corner)) then
    c(1) = corner(1)
  endif
  e(:) = 1
  if (present(edge_lengths)) then
    e(1) = edge_lengths(1)
  else
    e(1) = size(buf)
  endif
  if (present(unlim_dim_level)) then
    unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                             broadcast=bcast)
    if (unlim_dim_index .ne. 2) then
      call error("unlimited dimension must be the slowest varying dimension.")
    endif
    c(unlim_dim_index) = unlim_dim_level
  endif
  if (fileobj%is_root) then
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(buf)
      type is (integer(kind=int32))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (integer(kind=int64))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (real(kind=real32))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (real(kind=real64))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (character(len=*))
        ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.false.)
        if (ndims .ne. 2) then
          call error(trim(variable_name)//"must be 2d dimensional in netcdf file.")
        endif
        start(1) = 1
        start(2) = c(1)
        call get_variable_size(fileobj, variable_name, dimsizes, .false.)
        dimsizes(2) = dimsizes(2) - start(2) + 1
        call allocate_array(charbuf, dimsizes)
        charbuf(:,:) = ""
        err = nf90_get_var(fileobj%ncid, varid, charbuf, start=start, count=dimsizes)
        if (len(buf(1)) .lt. dimsizes(1)) then
          call error("character buffer is too small; increase length.")
        endif
        if (size(buf) .lt. dimsizes(2)) then
          call error("incorrect buffer array size.")
        endif
        do i = start(2), start(2)+dimsizes(2)-1
          sbuf = ""
          do j = 1, dimsizes(1)
            if (charbuf(j, i-start(2)+1) .eq. char(0)) then
              exit
            endif
            sbuf(j:j) = charbuf(j, i-start(2)+1)
          enddo
          call string_copy(buf(i-start(2)+1), sbuf)
        enddo
        deallocate(charbuf)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)
  endif
  if (bcast) then
    select type(buf)
      type is (integer(kind=int32))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (integer(kind=int64))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (real(kind=real32))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (real(kind=real64))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (character(len=*))
        call mpp_broadcast(buf, len(buf(1)), fileobj%io_root, pelist=fileobj%pelist)
      class default
        call error("unsupported type.")
    end select
  endif
end subroutine netcdf_read_data_1d


!> @brief Read in data from a variable in a netcdf file.
subroutine netcdf_read_data_2d(fileobj, variable_name, buf, unlim_dim_level, &
                               corner, edge_lengths, broadcast)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(inout) :: buf !< Array that the data
                                                 !! will be read into.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be read from.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be read
                                                              !! in each dimension.
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.

  logical :: bcast
  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer, dimension(3) :: c
  integer, dimension(3) :: e

  if (present(broadcast)) then
    bcast = broadcast
  else
    bcast = .true.
  endif
  c(:) = 1
  if (present(corner)) then
    c(1:2) = corner(1:2)
  endif
  e(:) = 1
  if (present(edge_lengths)) then
    e(1:2) = edge_lengths(1:2)
  else
    e(1:2) = shape(buf)
  endif
  if (present(unlim_dim_level)) then
    unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                             broadcast=bcast)
    if (unlim_dim_index .ne. 3) then
      call error("unlimited dimension must be the slowest varying dimension.")
    endif
    c(unlim_dim_index) = unlim_dim_level
  endif
  if (fileobj%is_root) then
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(buf)
      type is (integer(kind=int32))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (integer(kind=int64))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (real(kind=real32))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (real(kind=real64))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)
  endif
  if (bcast) then
    select type(buf)
      type is (integer(kind=int32))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (integer(kind=int64))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (real(kind=real32))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (real(kind=real64))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      class default
        call error("unsupported type.")
    end select
  endif
end subroutine netcdf_read_data_2d


!> @brief Read in data from a variable in a netcdf file.
subroutine netcdf_read_data_3d(fileobj, variable_name, buf, unlim_dim_level, &
                               corner, edge_lengths, broadcast)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(inout) :: buf !< Array that the data
                                                   !! will be read into.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be read from.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be read
                                                              !! in each dimension.
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.

  logical :: bcast
  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer, dimension(4) :: c
  integer, dimension(4) :: e

  if (present(broadcast)) then
    bcast = broadcast
  else
    bcast = .true.
  endif
  c(:) = 1
  if (present(corner)) then
    c(1:3) = corner(1:3)
  endif
  e(:) = 1
  if (present(edge_lengths)) then
    e(1:3) = edge_lengths(1:3)
  else
    e(1:3) = shape(buf)
  endif
  if (present(unlim_dim_level)) then
    unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                             broadcast=bcast)
    if (unlim_dim_index .ne. 4) then
      call error("unlimited dimension must be the slowest varying dimension.")
    endif
    c(unlim_dim_index) = unlim_dim_level
  endif
  if (fileobj%is_root) then
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(buf)
      type is (integer(kind=int32))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (integer(kind=int64))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (real(kind=real32))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (real(kind=real64))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)
  endif
  if (bcast) then
    select type(buf)
      type is (integer(kind=int32))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (integer(kind=int64))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (real(kind=real32))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (real(kind=real64))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      class default
        call error("unsupported type.")
    end select
  endif
end subroutine netcdf_read_data_3d


!> @brief Read in data from a variable in a netcdf file.
subroutine netcdf_read_data_4d(fileobj, variable_name, buf, unlim_dim_level, &
                               corner, edge_lengths, broadcast)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(inout) :: buf !< Array that the data
                                                     !! will be read into.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be read from.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be read
                                                              !! in each dimension.
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.

  logical :: bcast
  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer, dimension(5) :: c
  integer, dimension(5) :: e

  if (present(broadcast)) then
    bcast = broadcast
  else
    bcast = .true.
  endif
  c(:) = 1
  if (present(corner)) then
    c(1:4) = corner(1:4)
  endif
  e(:) = 1
  if (present(edge_lengths)) then
    e(1:4) = edge_lengths(1:4)
  else
    e(1:4) = shape(buf)
  endif
  if (present(unlim_dim_level)) then
    unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                             broadcast=bcast)
    if (unlim_dim_index .ne. 5) then
      call error("unlimited dimension must be the slowest varying dimension.")
    endif
    c(unlim_dim_index) = unlim_dim_level
  endif
  if (fileobj%is_root) then
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(buf)
      type is (integer(kind=int32))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (integer(kind=int64))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (real(kind=real32))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (real(kind=real64))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)
  endif
  if (bcast) then
    select type(buf)
      type is (integer(kind=int32))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (integer(kind=int64))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (real(kind=real32))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (real(kind=real64))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      class default
        call error("unsupported type.")
    end select
  endif
end subroutine netcdf_read_data_4d


!> @brief Read in data from a variable in a netcdf file.
subroutine netcdf_read_data_5d(fileobj, variable_name, buf, unlim_dim_level, &
                               corner, edge_lengths, broadcast)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(inout) :: buf !< Array that the data
                                                       !! will be read into.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(:), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be read from.
  integer, dimension(:), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be read
                                                              !! in each dimension.
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.

  logical :: bcast
  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer, dimension(6) :: c
  integer, dimension(6) :: e

  if (present(broadcast)) then
    bcast = broadcast
  else
    bcast = .true.
  endif
  c(:) = 1
  if (present(corner)) then
    c(1:5) = corner(1:5)
  endif
  e(:) = 1
  if (present(edge_lengths)) then
    e(1:5) = edge_lengths(1:5)
  else
    e(1:5) = shape(buf)
  endif
  if (present(unlim_dim_level)) then
    unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                             broadcast=bcast)
    if (unlim_dim_index .ne. 6) then
      call error("unlimited dimension must be the slowest varying dimension.")
    endif
    c(unlim_dim_index) = unlim_dim_level
  endif
  if (fileobj%is_root) then
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(buf)
      type is (integer(kind=int32))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (integer(kind=int64))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (real(kind=real32))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      type is (real(kind=real64))
        err = nf90_get_var(fileobj%ncid, varid, buf, start=c, count=e)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)
  endif
  if (bcast) then
    select type(buf)
      type is (integer(kind=int32))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (integer(kind=int64))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (real(kind=real32))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      type is (real(kind=real64))
        call mpp_broadcast(buf, size(buf), fileobj%io_root, pelist=fileobj%pelist)
      class default
        call error("unsupported type.")
    end select
  endif
end subroutine netcdf_read_data_5d
