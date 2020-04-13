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

!> @briefs Calculates a variable's checksum across all ranks in the current pelist.
!! @return A hex string containing the checksum.
function compute_global_checksum_2d(fileobj, variable_name, variable_data, is_decomposed) &
  result(chksum)

  class(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in) :: variable_data !< Data to be checksummed.
  logical, intent(out) :: is_decomposed !< Flag telling if the variable is decomposed.
  character(len=32) :: chksum

  integer :: xdim
  integer :: ydim
  type(domain2d), pointer :: io_domain
  integer :: xpos
  integer :: ypos
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  logical :: buffer_includes_halos
  logical :: extra_x
  logical :: extra_y
  integer, dimension(2) :: c
  integer, dimension(2) :: e
  integer(kind=int32), dimension(:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:), allocatable :: buf_real64
  integer(kind=int32) :: fill_int32
  integer(kind=int64) :: fill_int64
  real(kind=real32) :: fill_real32
  real(kind=real64) :: fill_real64
  integer(kind=int64) :: chksum_val

  is_decomposed = is_variable_domain_decomposed(fileobj, variable_name, &
                                                xindex=xdim, yindex=ydim, &
                                                xpos=xpos, ypos=ypos)
  if (.not. is_decomposed) then
    chksum = ""
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  call domain_offsets(size(variable_data, xdim), size(variable_data, ydim), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, &
                      yc_size, buffer_includes_halos, extra_x, extra_y)
  c(:) = 1
  if (buffer_includes_halos) then
    !Adjust if the input buffer has room for halos.
    c(xdim) = isc - isd + 1
    c(ydim) = jsc - jsd + 1
  endif
  e(:) = shape(variable_data)
  e(xdim) = xc_size
  e(ydim) = yc_size

  if (extra_x) then
    !Adjust sizes since compute domains overlap whe there are non-centered
    !domain position.
    e(xdim) = e(xdim) - 1
  endif

  if (extra_y) then
    !Adjust sizes since compute domains overlap whe there are non-centered
    !domain position.
    e(ydim) = e(ydim) - 1
  endif

  select type (variable_data)
    type is (integer(kind=int32))
      call allocate_array(buf_int32, e)
      call get_array_section(buf_int32, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_int32)) then
        chksum_val = mpp_chksum(buf_int32, mask_val=fill_int32)
      else
        chksum_val = mpp_chksum(buf_int32)
      endif
      deallocate(buf_int32)
    type is (integer(kind=int64))
      call allocate_array(buf_int64, e)
      call get_array_section(buf_int64, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_int64)) then
        chksum_val = mpp_chksum(buf_int64, mask_val=fill_int64)
      else
        chksum_val = mpp_chksum(buf_int64)
      endif
      deallocate(buf_int64)
   type is (real(kind=real32))
#ifdef OVERLOAD_R4
     call allocate_array(buf_real32, e)
     call get_array_section(buf_real32, variable_data, c, e)
     if (get_fill_value(fileobj, variable_name, fill_real32)) then
       chksum_val = mpp_chksum(buf_real32, mask_val=fill_real32)
     else
       chksum_val = mpp_chksum(buf_real32)
     endif
     deallocate(buf_real32)
#else
     call error("Compute_global_checksum: you are trying to use a real*4 without defining OVERLOAD_R4")
#endif
    type is (real(kind=real64))
      call allocate_array(buf_real64, e)
      call get_array_section(buf_real64, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_real64)) then
        chksum_val = mpp_chksum(buf_real64, mask_val=fill_real64)
      else
        chksum_val = mpp_chksum(buf_real64)
      endif
      deallocate(buf_real64)
    class default
      call error("Compute_global_checksum: unsupported type.")
  end select
  chksum = ""
  write(chksum, "(Z16)") chksum_val
end function compute_global_checksum_2d


!> @briefs Calculates a variable's checksum across all ranks in the current pelist.
!! @return A hex string containing the checksum.
function compute_global_checksum_3d(fileobj, variable_name, variable_data, is_decomposed) &
  result(chksum)

  class(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in) :: variable_data !< Data to be checksummed.
  logical, intent(out) :: is_decomposed !< Flag telling if the variable is decomposed.
  character(len=32) :: chksum

  integer :: xdim
  integer :: ydim
  type(domain2d), pointer :: io_domain
  integer :: xpos
  integer :: ypos
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  logical :: buffer_includes_halos
  logical :: extra_x
  logical :: extra_y
  integer, dimension(3) :: c
  integer, dimension(3) :: e
  integer(kind=int32), dimension(:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:), allocatable :: buf_real64
  integer(kind=int32) :: fill_int32
  integer(kind=int64) :: fill_int64
  real(kind=real32) :: fill_real32
  real(kind=real64) :: fill_real64
  integer(kind=int64) :: chksum_val

  is_decomposed = is_variable_domain_decomposed(fileobj, variable_name, &
                                                xindex=xdim, yindex=ydim, &
                                                xpos=xpos, ypos=ypos)
  if (.not. is_decomposed) then
    chksum = ""
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  call domain_offsets(size(variable_data, xdim), size(variable_data, ydim), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, &
                      yc_size, buffer_includes_halos, extra_x, extra_y)
  c(:) = 1
  if (buffer_includes_halos) then
    !Adjust if the input buffer has room for halos.
    c(xdim) = isc - isd + 1
    c(ydim) = jsc - jsd + 1
  endif
  e(:) = shape(variable_data)
  e(xdim) = xc_size
  e(ydim) = yc_size

  if (extra_x) then
    !Adjust sizes since compute domains overlap whe there are non-centered
    !domain position.
    e(xdim) = e(xdim) - 1
  endif

  if (extra_y) then
    !Adjust sizes since compute domains overlap whe there are non-centered
    !domain position.
    e(ydim) = e(ydim) - 1
  endif
  select type (variable_data)
    type is (integer(kind=int32))
      call allocate_array(buf_int32, e)
      call get_array_section(buf_int32, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_int32)) then
        chksum_val = mpp_chksum(buf_int32, mask_val=fill_int32)
      else
        chksum_val = mpp_chksum(buf_int32)
      endif
      deallocate(buf_int32)
    type is (integer(kind=int64))
      call allocate_array(buf_int64, e)
      call get_array_section(buf_int64, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_int64)) then
        chksum_val = mpp_chksum(buf_int64, mask_val=fill_int64)
      else
        chksum_val = mpp_chksum(buf_int64)
      endif
      deallocate(buf_int64)
   type is (real(kind=real32))
#ifdef OVERLOAD_R4
     call allocate_array(buf_real32, e)
     call get_array_section(buf_real32, variable_data, c, e)
     if (get_fill_value(fileobj, variable_name, fill_real32)) then
       chksum_val = mpp_chksum(buf_real32, mask_val=fill_real32)
     else
       chksum_val = mpp_chksum(buf_real32)
     endif
     deallocate(buf_real32)
#else
     call error("Compute_global_checksum: you are trying to use a real*4 without defining OVERLOAD_R4")
#endif
    type is (real(kind=real64))
      call allocate_array(buf_real64, e)
      call get_array_section(buf_real64, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_real64)) then
        chksum_val = mpp_chksum(buf_real64, mask_val=fill_real64)
      else
        chksum_val = mpp_chksum(buf_real64)
      endif
      deallocate(buf_real64)
    class default
      call error("Compute_global_checksum: unsupported type.")
  end select
  chksum = ""
  write(chksum, "(Z16)") chksum_val
end function compute_global_checksum_3d


!> @briefs Calculates a variable's checksum across all ranks in the current pelist.
!! @return A hex string containing the checksum.
function compute_global_checksum_4d(fileobj, variable_name, variable_data, is_decomposed) &
  result(chksum)

  class(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(in) :: variable_data !< Data to be checksummed.
  logical, intent(out) :: is_decomposed !< Flag telling if the variable is decomposed.
  character(len=32) :: chksum

  integer :: xdim
  integer :: ydim
  type(domain2d), pointer :: io_domain
  integer :: xpos
  integer :: ypos
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  logical :: buffer_includes_halos
  logical :: extra_x
  logical :: extra_y
  integer, dimension(4) :: c
  integer, dimension(4) :: e
  integer(kind=int32), dimension(:,:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:,:), allocatable :: buf_real64
  integer(kind=int32) :: fill_int32
  integer(kind=int64) :: fill_int64
  real(kind=real32) :: fill_real32
  real(kind=real64) :: fill_real64
  integer(kind=int64) :: chksum_val

  is_decomposed = is_variable_domain_decomposed(fileobj, variable_name, &
                                                xindex=xdim, yindex=ydim, &
                                                xpos=xpos, ypos=ypos)
  if (.not. is_decomposed) then
    chksum = ""
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  call domain_offsets(size(variable_data, xdim), size(variable_data, ydim), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, &
                      yc_size, buffer_includes_halos, extra_x, extra_y)
  c(:) = 1
  if (buffer_includes_halos) then
    !Adjust if the input buffer has room for halos.
    c(xdim) = isc - isd + 1
    c(ydim) = jsc - jsd + 1
  endif
  e(:) = shape(variable_data)
  e(xdim) = xc_size
  e(ydim) = yc_size

  if (extra_x) then
    !Adjust sizes since compute domains overlap whe there are non-centered
    !domain position.
    e(xdim) = e(xdim) - 1
  endif

  if (extra_y) then
    !Adjust sizes since compute domains overlap whe there are non-centered
    !domain position.
    e(ydim) = e(ydim) - 1
  endif

  select type (variable_data)
    type is (integer(kind=int32))
      call allocate_array(buf_int32, e)
      call get_array_section(buf_int32, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_int32)) then
        chksum_val = mpp_chksum(buf_int32, mask_val=fill_int32)
      else
        chksum_val = mpp_chksum(buf_int32)
      endif
      deallocate(buf_int32)
    type is (integer(kind=int64))
      call allocate_array(buf_int64, e)
      call get_array_section(buf_int64, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_int64)) then
        chksum_val = mpp_chksum(buf_int64, mask_val=fill_int64)
      else
        chksum_val = mpp_chksum(buf_int64)
      endif
      deallocate(buf_int64)
   type is (real(kind=real32))
#ifdef OVERLOAD_R4
     call allocate_array(buf_real32, e)
     call get_array_section(buf_real32, variable_data, c, e)
     if (get_fill_value(fileobj, variable_name, fill_real32)) then
       chksum_val = mpp_chksum(buf_real32, mask_val=fill_real32)
     else
       chksum_val = mpp_chksum(buf_real32)
     endif
     deallocate(buf_real32)
#else
     call error("Compute_global_checksum: you are trying to use a real*4 without defining OVERLOAD_R4")
#endif
    type is (real(kind=real64))
      call allocate_array(buf_real64, e)
      call get_array_section(buf_real64, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_real64)) then
        chksum_val = mpp_chksum(buf_real64, mask_val=fill_real64)
      else
        chksum_val = mpp_chksum(buf_real64)
      endif
      deallocate(buf_real64)
    class default
      call error("Compute_global_checksum: unsupported type.")
  end select
  chksum = ""
  write(chksum, "(Z16)") chksum_val
end function compute_global_checksum_4d
