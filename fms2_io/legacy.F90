!> @brief Recreate some legacy behavior.
module legacy_mod

use, intrinsic :: iso_fortran_env
use netcdf
use fms2_io_mod
use fms_io_utils_mod

implicit none
private


public :: axis_edges


contains


subroutine axis_edges(fileobj, name, edge_data)

  class(FmsNetcdfFile_t), intent(in) :: fileobj
  character(len=*), intent(in) :: name
  class(*), dimension(:), intent(out) :: edge_data

  integer :: ndims
  character(len=nf90_max_name) :: buffer
  integer, dimension(:), allocatable :: dim_sizes
  real(kind=real32), dimension(:), allocatable :: r32
  real(kind=real32), dimension(:,:), allocatable :: r322d
  real(kind=real64), dimension(:), allocatable :: r64
  real(kind=real64), dimension(:,:), allocatable :: r642d
  integer :: i
  integer :: n

  ndims = get_variable_num_dimensions(fileobj, name)
  allocate(dim_sizes(ndims))
  call get_variable_size(fileobj, name, dim_sizes)
  n = dim_sizes(1)
  if (size(edge_data) .ne. n+1) then
    call error("incorrect size of edge_data array.")
  endif
  deallocate(dim_sizes)

  buffer = ""
  if (variable_att_exists(fileobj, name, "edges")) then
    call get_variable_attribute(fileobj, name, "edges", buffer)
  elseif (variable_att_exists(fileobj, name, "bounds")) then
    call get_variable_attribute(fileobj, name, "bounds", buffer)
  endif
  if (trim(buffer) .ne. "") then
    ndims = get_variable_num_dimensions(fileobj, buffer)
    allocate(dim_sizes(ndims))
    call get_variable_size(fileobj, buffer, dim_sizes)
    if (size(dim_sizes) .eq. 1) then
      if (dim_sizes(1) .ne. n+1) then
        call error("incorrect size of edge data.")
      endif
      select type (edge_data)
        type is (real(kind=real32))
          call read_data(fileobj, buffer, edge_data)
        type is (real(kind=real64))
          call read_data(fileobj, buffer, edge_data)
        class default
          call error("unsupported kind.")
      end select
    elseif (size(dim_sizes) .eq. 2) then
      if (dim_sizes(1) .ne. 2) then
        call error("first dimension of edge must be of size 2")
      endif
      if (dim_sizes(2) .ne. n) then
        call error("incorrect size of edge data.")
      endif
      select type (edge_data)
        type is (real(kind=real32))
          allocate(r322d(dim_sizes(1), dim_sizes(2)))
          call read_data(fileobj, buffer, r322d)
          edge_data(1:dim_sizes(2)) = r322d(1,:)
          edge_data(dim_sizes(2)+1) = r322d(2,dim_sizes(2))
          deallocate(r322d)
        type is (real(kind=real64))
          allocate(r642d(dim_sizes(1), dim_sizes(2)))
          call read_data(fileobj, buffer, r642d)
          edge_data(1:dim_sizes(2)) = r642d(1,:)
          edge_data(dim_sizes(2)+1) = r642d(2,dim_sizes(2))
          deallocate(r642d)
        class default
          call error("unsupported kind.")
      end select
    endif
    deallocate(dim_sizes)
  else
    select type (edge_data)
      type is (real(kind=real32))
        allocate(r32(n))
        call read_data(fileobj, name, r32)
        do i = 2, n
          edge_data(i) = r32(i-1) + 0.5_real32*(r32(i) - r32(i-1))
        enddo
        edge_data(1) = r32(1) - 0.5_real32*(r32(2) - r32(1))
        if (abs(edge_data(1)) .lt. 1.e-10) then
          edge_data(1) = 0._real32
        endif
        edge_data(n+1) = r32(n) + 0.5_real32*(r32(n) - r32(n-1))
        deallocate(r32)
      type is (real(kind=real64))
        allocate(r64(n))
        call read_data(fileobj, name, r64)
        do i = 2, n
          edge_data(i) = r64(i-1) + 0.5_real64*(r64(i) - r64(i-1))
        enddo
        edge_data(1) = r64(1) - 0.5_real64*(r64(2) - r64(1))
        if (abs(edge_data(1)) .lt. 1.d-10) then
          edge_data(1) = 0._real64
        endif
        edge_data(n+1) = r64(n) + 0.5_real64*(r64(n) - r64(n-1))
        deallocate(r64)
      class default
        call error("unsupported kind.")
    end select
  endif
end subroutine axis_edges


end module legacy_mod
