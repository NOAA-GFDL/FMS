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
!> @ingroup fms_netcdf_domain_io_mod

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
  integer(kind=i4_kind), dimension(:,:), allocatable :: buf_i4_kind
  integer(kind=i8_kind), dimension(:,:), allocatable :: buf_i8_kind
  real(kind=r4_kind), dimension(:,:), allocatable :: buf_r4_kind
  real(kind=r8_kind), dimension(:,:), allocatable :: buf_r8_kind
  integer(kind=i4_kind) :: fill_i4_kind
  integer(kind=i8_kind) :: fill_i8_kind
  real(kind=r4_kind) :: fill_r4_kind
  real(kind=r8_kind) :: fill_r8_kind
  integer(kind=i8_kind) :: chksum_val

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
                      yc_size, buffer_includes_halos, extra_x, extra_y, &
                      msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))
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
    type is (integer(kind=i4_kind))
      call allocate_array(buf_i4_kind, e)
      call get_array_section(buf_i4_kind, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_i4_kind)) then
        chksum_val = mpp_chksum(buf_i4_kind, mask_val=fill_i4_kind)
      else
        chksum_val = mpp_chksum(buf_i4_kind)
      endif
      deallocate(buf_i4_kind)
    type is (integer(kind=i8_kind))
      call allocate_array(buf_i8_kind, e)
      call get_array_section(buf_i8_kind, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_i8_kind)) then
        chksum_val = mpp_chksum(buf_i8_kind, mask_val=fill_i8_kind)
      else
        chksum_val = mpp_chksum(buf_i8_kind)
      endif
      deallocate(buf_i8_kind)
   type is (real(kind=r4_kind))
     call allocate_array(buf_r4_kind, e)
     call get_array_section(buf_r4_kind, variable_data, c, e)
     if (get_fill_value(fileobj, variable_name, fill_r4_kind)) then
       chksum_val = mpp_chksum(buf_r4_kind, mask_val=fill_r4_kind)
     else
       chksum_val = mpp_chksum(buf_r4_kind)
     endif
     deallocate(buf_r4_kind)
   type is (real(kind=r8_kind))
      call allocate_array(buf_r8_kind, e)
      call get_array_section(buf_r8_kind, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_r8_kind)) then
        chksum_val = mpp_chksum(buf_r8_kind, mask_val=fill_r8_kind)
      else
        chksum_val = mpp_chksum(buf_r8_kind)
      endif
      deallocate(buf_r8_kind)
    class default
      call error("unsupported variable type: compute_global_checksum_2d: file: "//trim(fileobj%path)//" variable:"//trim(variable_name))
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
  integer(kind=i4_kind), dimension(:,:,:), allocatable :: buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:), allocatable :: buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:), allocatable :: buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:), allocatable :: buf_r8_kind
  integer(kind=i4_kind) :: fill_i4_kind
  integer(kind=i8_kind) :: fill_i8_kind
  real(kind=r4_kind) :: fill_r4_kind
  real(kind=r8_kind) :: fill_r8_kind
  integer(kind=i8_kind) :: chksum_val

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
                      yc_size, buffer_includes_halos, extra_x, extra_y, &
                      msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))
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
    type is (integer(kind=i4_kind))
      call allocate_array(buf_i4_kind, e)
      call get_array_section(buf_i4_kind, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_i4_kind)) then
        chksum_val = mpp_chksum(buf_i4_kind, mask_val=fill_i4_kind)
      else
        chksum_val = mpp_chksum(buf_i4_kind)
      endif
      deallocate(buf_i4_kind)
    type is (integer(kind=i8_kind))
      call allocate_array(buf_i8_kind, e)
      call get_array_section(buf_i8_kind, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_i8_kind)) then
        chksum_val = mpp_chksum(buf_i8_kind, mask_val=fill_i8_kind)
      else
        chksum_val = mpp_chksum(buf_i8_kind)
      endif
      deallocate(buf_i8_kind)
   type is (real(kind=r4_kind))
     call allocate_array(buf_r4_kind, e)
     call get_array_section(buf_r4_kind, variable_data, c, e)
     if (get_fill_value(fileobj, variable_name, fill_r4_kind)) then
       chksum_val = mpp_chksum(buf_r4_kind, mask_val=fill_r4_kind)
     else
       chksum_val = mpp_chksum(buf_r4_kind)
     endif
     deallocate(buf_r4_kind)
   type is (real(kind=r8_kind))
      call allocate_array(buf_r8_kind, e)
      call get_array_section(buf_r8_kind, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_r8_kind)) then
        chksum_val = mpp_chksum(buf_r8_kind, mask_val=fill_r8_kind)
      else
        chksum_val = mpp_chksum(buf_r8_kind)
      endif
      deallocate(buf_r8_kind)
    class default
      call error("unsupported variable type: compute_global_checksum_3d: file: "//trim(fileobj%path)//" variable:"//trim(variable_name))
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
  integer(kind=i4_kind), dimension(:,:,:,:), allocatable :: buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:,:), allocatable :: buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:,:), allocatable :: buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:,:), allocatable :: buf_r8_kind
  integer(kind=i4_kind) :: fill_i4_kind
  integer(kind=i8_kind) :: fill_i8_kind
  real(kind=r4_kind) :: fill_r4_kind
  real(kind=r8_kind) :: fill_r8_kind
  integer(kind=i8_kind) :: chksum_val

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
                      yc_size, buffer_includes_halos, extra_x, extra_y, &
                      msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))
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
    type is (integer(kind=i4_kind))
      call allocate_array(buf_i4_kind, e)
      call get_array_section(buf_i4_kind, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_i4_kind)) then
        chksum_val = mpp_chksum(buf_i4_kind, mask_val=fill_i4_kind)
      else
        chksum_val = mpp_chksum(buf_i4_kind)
      endif
      deallocate(buf_i4_kind)
    type is (integer(kind=i8_kind))
      call allocate_array(buf_i8_kind, e)
      call get_array_section(buf_i8_kind, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_i8_kind)) then
        chksum_val = mpp_chksum(buf_i8_kind, mask_val=fill_i8_kind)
      else
        chksum_val = mpp_chksum(buf_i8_kind)
      endif
      deallocate(buf_i8_kind)
   type is (real(kind=r4_kind))
     call allocate_array(buf_r4_kind, e)
     call get_array_section(buf_r4_kind, variable_data, c, e)
     if (get_fill_value(fileobj, variable_name, fill_r4_kind)) then
       chksum_val = mpp_chksum(buf_r4_kind, mask_val=fill_r4_kind)
     else
       chksum_val = mpp_chksum(buf_r4_kind)
     endif
     deallocate(buf_r4_kind)
   type is (real(kind=r8_kind))
      call allocate_array(buf_r8_kind, e)
      call get_array_section(buf_r8_kind, variable_data, c, e)
      if (get_fill_value(fileobj, variable_name, fill_r8_kind)) then
        chksum_val = mpp_chksum(buf_r8_kind, mask_val=fill_r8_kind)
      else
        chksum_val = mpp_chksum(buf_r8_kind)
      endif
      deallocate(buf_r8_kind)
    class default
      call error("unsupported variable type: compute_global_checksum_4d: file: "//trim(fileobj%path)//" variable:"//trim(variable_name))
  end select
  chksum = ""
  write(chksum, "(Z16)") chksum_val
end function compute_global_checksum_4d
