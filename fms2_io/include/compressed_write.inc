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
!> @file
!> @brief Compressed write routines for @ref write_data interface

!> @addtogroup fms2_io_mod
!> @{

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

  class(FmsNetcdfFile_t), intent(in)           :: fileobj         !< File object.
  character(len=*),       intent(in)           :: variable_name   !< Variable name.
  class(*), dimension(:), intent(in)           :: cdata           !< Compressed data that
                                                                  !! will be gathered and
                                                                  !! written to the
                                                                  !! netcdf file.
  integer,                intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                                  !! dimension.
  integer,  dimension(1), intent(in), optional :: corner          !< Array of starting
                                                                  !! indices describing
                                                                  !! where the data
                                                                  !! will be written to.
  integer,  dimension(1), intent(in), optional :: edge_lengths    !< The number of
                                                                  !! elements that
                                                                  !! will be written
                                                                  !! in each dimension.

  integer, dimension(2) :: compressed_dim_index !< index of the compressed dimension
                                                !! compressed_dim_index(1) relative to cdata
                                                !! compressed_dim_index(2) relative to the fileobj
  integer, dimension(1) :: e                    !< "edges" number of points to read

  integer(kind=i4_kind), dimension(:), allocatable :: buf_i4_kind !< Global buffer of data
  integer(kind=i8_kind), dimension(:), allocatable :: buf_i8_kind !< Global buffer of data
  real   (kind=r4_kind), dimension(:), allocatable :: buf_r4_kind !< Global buffer of data
  real   (kind=r8_kind), dimension(:), allocatable :: buf_r8_kind !< Global buffer of data

  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message

  append_error_msg = "compressed_write_1d: file:"//trim(fileobj%path)//" variable:"//trim(variable_name)

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)

  if (compressed_dim_index(1) .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner, &
                           edge_lengths=edge_lengths)
    return
  endif

  e(:) = shape(cdata)
  !The root pe creates a buffer big enough to store the data:
  if (fileobj%is_root) then
    e(compressed_dim_index(1)) = sum(fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems)
    select type(cdata)
      type is (integer(kind=i4_kind))
        call allocate_array(buf_i4_kind, e)
      type is (integer(kind=i8_kind))
        call allocate_array(buf_i8_kind, e)
      type is (real(kind=r4_kind))
        call allocate_array(buf_r4_kind, e)
      type is (real(kind=r8_kind))
        call allocate_array(buf_r8_kind, e)
      class default
        call error("unsupported variable type: "//trim(append_error_msg))
     end select
  endif

  !Gather the data onto the I/O root and write it out.
  select type(cdata)
    type is (integer(kind=i4_kind))
      call mpp_gather(cdata, size(cdata), buf_i4_kind, fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems, &
        fileobj%pelist)
      call netcdf_write_data(fileobj, variable_name, buf_i4_kind, &
        unlim_dim_level=unlim_dim_level)
    type is (integer(kind=i8_kind))
      call mpp_gather(cdata, size(cdata), buf_i8_kind, fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems, &
        fileobj%pelist)
      call netcdf_write_data(fileobj, variable_name, buf_i8_kind, &
        unlim_dim_level=unlim_dim_level)
    type is (real(kind=r4_kind))
      call mpp_gather(cdata, size(cdata), buf_r4_kind, fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems, &
        fileobj%pelist)
      call netcdf_write_data(fileobj, variable_name, buf_r4_kind, &
        unlim_dim_level=unlim_dim_level)
    type is (real(kind=r8_kind))
      call mpp_gather(cdata, size(cdata), buf_r8_kind, fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems, &
        fileobj%pelist)
      call netcdf_write_data(fileobj, variable_name, buf_r8_kind, &
        unlim_dim_level=unlim_dim_level)
    class default
      call error("unsupported variable type: "//trim(append_error_msg))
  end select
  if (fileobj%is_root) then
    if (allocated(buf_i4_kind)) deallocate(buf_i4_kind)
    if (allocated(buf_i8_kind)) deallocate(buf_i8_kind)
    if (allocated(buf_r4_kind)) deallocate(buf_r4_kind)
    if (allocated(buf_r8_kind)) deallocate(buf_r8_kind)
  endif
end subroutine compressed_write_1d


!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_2d(fileobj, variable_name, cdata, unlim_dim_level, &
                               corner, edge_lengths)

  class(FmsNetcdfFile_t),   intent(in)           :: fileobj         !< File object.
  character(len=*),         intent(in)           :: variable_name   !< Variable name.
  class(*), dimension(:,:), intent(in)           :: cdata           !< Compressed data that
                                                                    !! will be gathered and
                                                                    !! written to the
                                                                    !! netcdf file.
  integer,                  intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                                    !! dimension.
  integer,  dimension(2),   intent(in), optional :: corner          !< Array of starting
                                                                    !! indices describing
                                                                    !! where the data
                                                                    !! will be written to.
  integer,  dimension(2),   intent(in), optional :: edge_lengths    !< The number of
                                                                    !! elements that
                                                                    !! will be written
                                                                    !! in each dimension.

  integer, dimension(2) :: compressed_dim_index !< index of the compressed dimension
                                                !! compressed_dim_index(1) relative to cdata
                                                !! compressed_dim_index(2) relative to the fileobj
  integer, dimension(2) :: c                    !! corners of the data to read
  integer, dimension(2) :: e                    !< "edges" number of points to read

  integer(kind=i4_kind), dimension(:,:), allocatable :: buf_i4_kind !< Global buffer of data
  integer(kind=i8_kind), dimension(:,:), allocatable :: buf_i8_kind !< Global buffer of data
  real   (kind=r4_kind), dimension(:,:), allocatable :: buf_r4_kind !< Global buffer of data
  real   (kind=r8_kind), dimension(:,:), allocatable :: buf_r8_kind !< Global buffer of data

  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message

  integer :: index(1) !< index of the PE in the pelist

  integer :: is !< Starting index of the first dimension
  integer :: js !< Ending index of the first dimension
  integer :: ie !< Starting index of the second dimension
  integer :: je !< Ending index of the second dimension

  append_error_msg = "compressed_write_2d: file:"//trim(fileobj%path)//" variable:"//trim(variable_name)

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)
  if (compressed_dim_index(1) .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner, &
                           edge_lengths=edge_lengths)
    return
  endif

  e(:) = shape(cdata)
  !The root pe creates a buffer big enough to store the data:
  if (fileobj%is_root) then
    e(compressed_dim_index(1)) = sum(fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems)
    select type(cdata)
      type is (integer(kind=i4_kind))
        call allocate_array(buf_i4_kind, e)
      type is (integer(kind=i8_kind))
        call allocate_array(buf_i8_kind, e)
      type is (real(kind=r4_kind))
        call allocate_array(buf_r4_kind, e)
      type is (real(kind=r8_kind))
        call allocate_array(buf_r8_kind, e)
      class default
        call error("unsupported variable type: "//trim(append_error_msg))
     end select
  endif

  c(:) = 1
  index = FINDLOC(fileobj%pelist, mpp_pe())

  c(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_corner(index(1))
  e(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems(index(1))

  if (compressed_dim_index(1) .eq. 1) then
    is = c(compressed_dim_index(1))
    ie = c(compressed_dim_index(1)) + e(compressed_dim_index(1)) - 1
    js = 1
    je = size(cdata,2)
  else
    js = c(compressed_dim_index(1))
    je = c(compressed_dim_index(1)) + e(compressed_dim_index(1)) - 1
    is = 1
    ie = size(cdata,2)
  endif

  select type(cdata)
    type is (integer(kind=i4_kind))
      call mpp_gather(is, ie, js, je, fileobj%pelist, cdata, buf_i4_kind, fileobj%is_root)
      call netcdf_write_data(fileobj, variable_name, buf_i4_kind, &
        unlim_dim_level=unlim_dim_level)
    type is (integer(kind=i8_kind))
      call mpp_gather(is, ie, js, je, fileobj%pelist, cdata, buf_i8_kind, fileobj%is_root)
      call netcdf_write_data(fileobj, variable_name, buf_i8_kind, &
        unlim_dim_level=unlim_dim_level)
    type is (real(kind=r4_kind))
      call mpp_gather(is, ie, js, je, fileobj%pelist, cdata, buf_r4_kind, fileobj%is_root)
      call netcdf_write_data(fileobj, variable_name, buf_r4_kind, &
        unlim_dim_level=unlim_dim_level)
    type is (real(kind=r8_kind))
      call mpp_gather(is, ie, js, je, fileobj%pelist, cdata, buf_r8_kind, fileobj%is_root)
      call netcdf_write_data(fileobj, variable_name, buf_r8_kind, &
        unlim_dim_level=unlim_dim_level)
    class default
      call error("unsupported variable type: "//trim(append_error_msg))
  end select

  if (fileobj%is_root) then
    if (allocated(buf_i4_kind)) deallocate(buf_i4_kind)
    if (allocated(buf_i8_kind)) deallocate(buf_i8_kind)
    if (allocated(buf_r4_kind)) deallocate(buf_r4_kind)
    if (allocated(buf_r8_kind)) deallocate(buf_r8_kind)
  endif
end subroutine compressed_write_2d


!> @brief For variables without a compressed dimension, this routine simply wraps
!!        netcdf_write data.  If the variable does have a compressed axis, the I/O root
!!        gathers the data from the rest of the ranks and then writes the combined data to
!!        the netcdf file.
subroutine compressed_write_3d(fileobj, variable_name, cdata, unlim_dim_level, &
                               corner, edge_lengths)

  class(FmsNetcdfFile_t),     intent(in)           :: fileobj         !< File object.
  character(len=*),           intent(in)           :: variable_name   !< Variable name.
  class(*), contiguous,       intent(in), target   :: cdata(:,:,:)    !< Compressed data that
                                                                      !! will be gathered and
                                                                      !! written to the
                                                                      !! netcdf file.
  integer,                    intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                                      !! dimension.
  integer,  dimension(3),     intent(in), optional :: corner          !< Array of starting
                                                                      !! indices describing
                                                                      !! where the data
                                                                      !! will be written to.
  integer,  dimension(3),     intent(in), optional :: edge_lengths    !< The number of
                                                                      !! elements that
                                                                      !! will be written
                                                                      !! in each dimension.

  integer, dimension(2) :: compressed_dim_index !< index of the compressed dimension
                                                !! compressed_dim_index(1) relative to cdata
                                                !! compressed_dim_index(2) relative to the fileobj
  integer, dimension(3) :: c                    !! corners of the data to read
  integer, dimension(3) :: e                    !< "edges" number of points to read

  integer(kind=i4_kind), dimension(:,:,:), allocatable :: buf_i4_kind !< Global buffer of data
  integer(kind=i8_kind), dimension(:,:,:), allocatable :: buf_i8_kind !< Global buffer of data
  real   (kind=r4_kind), dimension(:,:,:), allocatable :: buf_r4_kind !< Global buffer of data
  real   (kind=r8_kind), dimension(:,:,:), allocatable :: buf_r8_kind !< Global buffer of data

  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message
  class(*), pointer :: cdata_dummy(:,:,:,:)

  integer :: index(1) !< index of the PE in the pelist

  integer :: is !< Starting index of the first dimension
  integer :: js !< Ending index of the first dimension
  integer :: ie !< Starting index of the second dimension
  integer :: je !< Ending index of the second dimension

  append_error_msg = "compressed_write_3d: file:"//trim(fileobj%path)//" variable:"//trim(variable_name)

  compressed_dim_index = get_variable_compressed_dimension_index(fileobj, variable_name)
  if (compressed_dim_index(1) .eq. dimension_not_found) then
    call netcdf_write_data(fileobj, variable_name, cdata, &
                           unlim_dim_level=unlim_dim_level, corner=corner, &
                           edge_lengths=edge_lengths)
    return
  else if (compressed_dim_index(1) .eq. 3) then
    cdata_dummy(1:size(cdata,1), 1:size(cdata,2), 1:size(cdata,3), 1:1) => cdata(:,:,:)
    call compressed_write_4d(fileobj, variable_name, cdata_dummy, unlim_dim_level)
  endif

  e(:) = shape(cdata)
  !The root pe creates a buffer big enough to store the data:
  if (fileobj%is_root) then
    e(compressed_dim_index(1)) = sum(fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems)
    select type(cdata)
      type is (integer(kind=i4_kind))
        call allocate_array(buf_i4_kind, e)
      type is (integer(kind=i8_kind))
        call allocate_array(buf_i8_kind, e)
      type is (real(kind=r4_kind))
        call allocate_array(buf_r4_kind, e)
      type is (real(kind=r8_kind))
        call allocate_array(buf_r8_kind, e)
      class default
        call error("unsupported variable type: "//trim(append_error_msg))
     end select
  endif

  c(:) = 1
  index = FINDLOC(fileobj%pelist, mpp_pe())

  c(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_corner(index(1))
  e(compressed_dim_index(1)) = fileobj%compressed_dims(compressed_dim_index(2))%npes_nelems(index(1))

  if (compressed_dim_index(1) .eq. 1) then
    is = c(compressed_dim_index(1))
    ie = c(compressed_dim_index(1)) + e(compressed_dim_index(1)) - 1
    js = 1
    je = size(cdata,2)
  else
    js = c(compressed_dim_index(1))
    je = c(compressed_dim_index(1)) + e(compressed_dim_index(1)) - 1
    is = 1
    ie = size(cdata,2)
  endif

  select type(cdata)
    type is (integer(kind=i4_kind))
      call mpp_gather(is, ie, js, je, size(cdata,3), &
        fileobj%pelist, cdata, buf_i4_kind, fileobj%is_root)
      call netcdf_write_data(fileobj, variable_name, buf_i4_kind, &
        unlim_dim_level=unlim_dim_level)
    type is (integer(kind=i8_kind))
      call mpp_gather(is, ie, js, je, size(cdata,3), &
        fileobj%pelist, cdata, buf_i8_kind, fileobj%is_root)
      call netcdf_write_data(fileobj, variable_name, buf_i8_kind, &
        unlim_dim_level=unlim_dim_level)
    type is (real(kind=r4_kind))
      call mpp_gather(is, ie, js, je, size(cdata,3), &
        fileobj%pelist, cdata, buf_r4_kind, fileobj%is_root)
      call netcdf_write_data(fileobj, variable_name, buf_r4_kind, &
        unlim_dim_level=unlim_dim_level)
    type is (real(kind=r8_kind))
      call mpp_gather(is, ie, js, je, size(cdata,3), &
        fileobj%pelist, cdata, buf_r8_kind, fileobj%is_root)
      call netcdf_write_data(fileobj, variable_name, buf_r8_kind, &
        unlim_dim_level=unlim_dim_level)
    class default
      call error("unsupported variable type: "//trim(append_error_msg))
  end select

  if (fileobj%is_root) then
    if (allocated(buf_i4_kind)) deallocate(buf_i4_kind)
    if (allocated(buf_i8_kind)) deallocate(buf_i8_kind)
    if (allocated(buf_r4_kind)) deallocate(buf_r4_kind)
    if (allocated(buf_r8_kind)) deallocate(buf_r8_kind)
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
  integer(kind=i4_kind), dimension(:,:,:,:), allocatable :: buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:,:), allocatable :: buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:,:), allocatable :: buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:,:), allocatable :: buf_r8_kind
  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message

  append_error_msg = "compressed_write_4d: file:"//trim(fileobj%path)//" variable:"//trim(variable_name)

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
          type is (integer(kind=i4_kind))
            call allocate_array(buf_i4_kind, e)
            call mpp_recv(buf_i4_kind, size(buf_i4_kind), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_i4_kind, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_i4_kind)
          type is (integer(kind=i8_kind))
            call allocate_array(buf_i8_kind, e)
            call mpp_recv(buf_i8_kind, size(buf_i8_kind), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_i8_kind, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_i8_kind)
          type is (real(kind=r4_kind))
            call allocate_array(buf_r4_kind, e)
            call mpp_recv(buf_r4_kind, size(buf_r4_kind), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_r4_kind, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_r4_kind)
          type is (real(kind=r8_kind))
            call allocate_array(buf_r8_kind, e)
            call mpp_recv(buf_r8_kind, size(buf_r8_kind), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_r8_kind, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_r8_kind)
          class default
            call error("unsupported variable type: "//trim(append_error_msg))
        end select
      endif
    enddo
  else
    select type(cdata)
      type is (integer(kind=i4_kind))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (integer(kind=i8_kind))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=r4_kind))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=r8_kind))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      class default
        call error("unsupported variable type: "//trim(append_error_msg))
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
  integer(kind=i4_kind), dimension(:,:,:,:,:), allocatable :: buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:,:,:), allocatable :: buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:,:,:), allocatable :: buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:,:,:), allocatable :: buf_r8_kind
  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message

  append_error_msg = "compressed_write_5d: file:"//trim(fileobj%path)//" variable:"//trim(variable_name)

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
          type is (integer(kind=i4_kind))
            call allocate_array(buf_i4_kind, e)
            call mpp_recv(buf_i4_kind, size(buf_i4_kind), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_i4_kind, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_i4_kind)
          type is (integer(kind=i8_kind))
            call allocate_array(buf_i8_kind, e)
            call mpp_recv(buf_i8_kind, size(buf_i8_kind), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_i8_kind, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_i8_kind)
          type is (real(kind=r4_kind))
            call allocate_array(buf_r4_kind, e)
            call mpp_recv(buf_r4_kind, size(buf_r4_kind), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_r4_kind, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_r4_kind)
          type is (real(kind=r8_kind))
            call allocate_array(buf_r8_kind, e)
            call mpp_recv(buf_r8_kind, size(buf_r8_kind), fileobj%pelist(i), block=.true.)
            call netcdf_write_data(fileobj, variable_name, buf_r8_kind, &
                                   unlim_dim_level=unlim_dim_level, corner=c, &
                                   edge_lengths=e)
            deallocate(buf_r8_kind)
          class default
            call error("unsupported variable type: "//trim(append_error_msg))
        end select
      endif
    enddo
  else
    select type(cdata)
      type is (integer(kind=i4_kind))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (integer(kind=i8_kind))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=r4_kind))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      type is (real(kind=r8_kind))
        call mpp_send(cdata, size(cdata), fileobj%io_root)
      class default
        call error("unsupported variable type: "//trim(append_error_msg))
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
!> @}
