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

!> @brief Write data to a variable in a netcdf file.
subroutine netcdf_write_data_0d(fileobj, variable_name, variable_data, unlim_dim_level, &
                                corner)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*),  intent(in) :: variable_data !< Data that will be written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.
  integer, intent(in), optional :: corner !< Array of starting
                                          !! indices describing
                                          !! where the data
                                          !! will be written to.

  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer, dimension(1) :: c
  integer :: ndims
  integer, dimension(:), allocatable :: start
  integer, dimension(:), allocatable :: dimsizes
  character, dimension(:), allocatable :: charbuf
  integer :: i
  integer :: tlen

  if (fileobj%is_root) then
    c(:) = 1
    if (present(corner)) then
      c(1) = corner
    endif
    if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                               broadcast=.false.)
      if (unlim_dim_index .ne. 1) then
        call error("unlimited dimension must be the slowest varying dimension.")
      endif
      c(unlim_dim_index) = unlim_dim_level
    endif
    call set_netcdf_mode(fileobj%ncid, data_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(variable_data)
      type is (integer(kind=int32))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c)
      type is (integer(kind=int64))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c)
      type is (real(kind=real32))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c)
      type is (real(kind=real64))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c)
      type is (character(len=*))
        ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.false.)
        if (ndims .ne. 1) then
          call error("currently only scalar and 1d character writes are supported.")
        endif
        allocate(start(ndims))
        start(:) = 1
        allocate(dimsizes(ndims))
        call get_variable_size(fileobj, variable_name, dimsizes, .false.)
        call allocate_array(charbuf, dimsizes)
        charbuf(:) = ""
        tlen = len_trim(variable_data)
        if (tlen .gt. dimsizes(1)) then
          call error("character buffer is too big; decrease length.")
        endif
        do i = 1, tlen
          charbuf(i) = variable_data(i:i)
        enddo
        err = nf90_put_var(fileobj%ncid, varid, charbuf, start=start, count=dimsizes)
        deallocate(charbuf)
        deallocate(dimsizes)
        deallocate(start)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)
  endif
end subroutine netcdf_write_data_0d


!> @brief Write data to a variable in a netcdf file.
subroutine netcdf_write_data_1d(fileobj, variable_name, variable_data, unlim_dim_level, &
                                corner, edge_lengths)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(in) :: variable_data !< Data that will be written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.
  integer, dimension(1), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(1), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer, dimension(2) :: c
  integer, dimension(2) :: e
  integer :: ndims
  integer, dimension(:), allocatable :: start
  integer, dimension(:), allocatable :: dimsizes
  character, dimension(:,:), allocatable :: charbuf
  character(len=1024) :: sbuf
  integer :: i
  integer :: j
  integer :: tlen

  if (fileobj%is_root) then
    c(:) = 1
    if (present(corner)) then
      c(1:1) = corner(:)
    endif
    e(:) = 1
    if (present(edge_lengths)) then
      e(1:1) = edge_lengths(:)
    else
      e(1:1) = shape(variable_data)
    endif
    if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                               broadcast=.false.)
      if (unlim_dim_index .ne. 2) then
        call error("unlimited dimension must be the slowest varying dimension.")
      endif
      c(unlim_dim_index) = unlim_dim_level
    endif
    call set_netcdf_mode(fileobj%ncid, data_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(variable_data)
      type is (integer(kind=int32))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (integer(kind=int64))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=real32))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=real64))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (character(len=*))
        ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.false.)
        if (ndims .ne. 2) then
          call error("currently only scalar and 1d character writes are supported.")
        endif
        allocate(start(ndims))
        start(:) = 1
        allocate(dimsizes(ndims))
        call get_variable_size(fileobj, variable_name, dimsizes, .false.)
        call allocate_array(charbuf, dimsizes)
        charbuf(:,:) = ""
        tlen = len(variable_data(1))
        if (tlen .gt. dimsizes(1)) then
          call error("character buffer is too big; decrease length.")
        endif
        if (size(variable_data) .ne. dimsizes(2)) then
          call error("incorrect size of variable_data array.")
        endif
        do j = 1, dimsizes(2)
          call string_copy(sbuf, variable_data(j))
          do i = 1, tlen
            charbuf(i,j) = sbuf(i:i)
          enddo
        enddo
        err = nf90_put_var(fileobj%ncid, varid, charbuf, start=start, count=dimsizes)
        deallocate(charbuf)
        deallocate(dimsizes)
        deallocate(start)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)
  endif
end subroutine netcdf_write_data_1d


!> @brief Write data to a variable in a netcdf file.
subroutine netcdf_write_data_2d(fileobj, variable_name, variable_data, unlim_dim_level, &
                                corner, edge_lengths)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in) :: variable_data !< Data that will be written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.
  integer, dimension(2), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(2), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer,dimension(3) :: c
  integer, dimension(3) :: e

  if (fileobj%is_root) then
    c(:) = 1
    if (present(corner)) then
      c(1:2) = corner(:)
    endif
    e(:) = 1
    if (present(edge_lengths)) then
      e(1:2) = edge_lengths(:)
    else
      e(1:2) = shape(variable_data)
    endif
    if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                               broadcast=.false.)
      if (unlim_dim_index .ne. 3) then
        call error("unlimited dimension must be the slowest varying dimension.")
      endif
      c(unlim_dim_index) = unlim_dim_level
    endif
    call set_netcdf_mode(fileobj%ncid, data_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(variable_data)
      type is (integer(kind=int32))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (integer(kind=int64))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=real32))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=real64))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c ,count=e)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)
  endif
end subroutine netcdf_write_data_2d


!> @brief Write data to a variable in a netcdf file.
subroutine netcdf_write_data_3d(fileobj, variable_name, variable_data, unlim_dim_level, &
                                corner, edge_lengths)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in) :: variable_data !< Data that will be written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(3), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(3), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer, dimension(4) :: c
  integer, dimension(4) :: e

  if (fileobj%is_root) then
    c(:) = 1
    if (present(corner)) then
      c(1:3) = corner(:)
    endif
    e(:) = 1
    if (present(edge_lengths)) then
      e(1:3) = edge_lengths(:)
    else
      e(1:3) = shape(variable_data)
    endif
    if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                               broadcast=.false.)
      if (unlim_dim_index .ne. 4) then
        call error("unlimited dimension must be the slowest varying dimension.")
      endif
      c(unlim_dim_index) = unlim_dim_level
    endif
    call set_netcdf_mode(fileobj%ncid, data_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(variable_data)
      type is (integer(kind=int32))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (integer(kind=int64))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=real32))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=real64))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)
  endif
end subroutine netcdf_write_data_3d


!> @brief Write data to a variable in a netcdf file.
subroutine netcdf_write_data_4d(fileobj, variable_name, variable_data, unlim_dim_level, &
                                corner, edge_lengths)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(in) :: variable_data !< Data that will be written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.
  integer, dimension(4), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(4), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer,dimension(5) :: c
  integer, dimension(5) :: e

  if (fileobj%is_root) then
    c(:) = 1
    if (present(corner)) then
      c(1:4) = corner(:)
    endif
    e(:) = 1
    if (present(edge_lengths)) then
      e(1:4) = edge_lengths(:)
    else
      e(1:4) = shape(variable_data)
    endif
    if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                               broadcast=.false.)
     if  (unlim_dim_index /= -1) then
       c(unlim_dim_index) = unlim_dim_level
     endif
    endif
    call set_netcdf_mode(fileobj%ncid, data_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(variable_data)
      type is (integer(kind=int32))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (integer(kind=int64))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=real32))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=real64))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)
  endif
end subroutine netcdf_write_data_4d


!> @brief Write data to a variable in a netcdf file.
subroutine netcdf_write_data_5d(fileobj, variable_name, variable_data, unlim_dim_level, &
                                corner, edge_lengths)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(in), target :: variable_data !< Data that will be
                                                                      !! written.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  integer, dimension(5), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(5), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer :: err
  integer :: varid
  integer :: unlim_dim_index
  integer,dimension(6) :: c
  integer, dimension(6) :: e

  if (fileobj%is_root) then
    c(:) = 1
    if (present(corner)) then
      c(1:5) = corner(:)
    endif
    e(:) = 1
    if (present(edge_lengths)) then
      e(1:5) = edge_lengths(:)
    else
      e(1:5) = shape(variable_data)
    endif
    if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, &
                                                               broadcast=.false.)
      if (unlim_dim_index .ne. 6) then
        call error("unlimited dimension must be the slowest varying dimension.")
      endif
      c(unlim_dim_index) = unlim_dim_level
    endif
    call set_netcdf_mode(fileobj%ncid, data_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(variable_data)
      type is (integer(kind=int32))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (integer(kind=int64))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=real32))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      type is (real(kind=real64))
        err = nf90_put_var(fileobj%ncid, varid, variable_data, start=c, count=e)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)
  endif
end subroutine netcdf_write_data_5d
