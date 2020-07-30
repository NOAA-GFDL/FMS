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

!> @brief Domain-specific I/O wrappers.
module fms_netcdf_domain_io_mod
use, intrinsic :: iso_fortran_env
use netcdf
use mpp_mod
use mpp_domains_mod
use fms_io_utils_mod
use netcdf_io_mod
implicit none
private


!Module constants.
integer, parameter :: no_domain_decomposed_dimension = 0
integer, parameter, public :: max_num_domain_decomposed_dims = 10
integer, parameter :: variable_not_found = 0
integer, parameter :: default_domain_position = center
character(len=16), parameter :: domain_pos_att = "domain_position"
character(len=16), parameter :: domain_axis_att_name = "domain_axis"
character(len=16), parameter :: x = "x"
character(len=16), parameter :: y = "y"


!> @brief Domain variable.
type :: DomainDimension_t
  character(len=nf90_max_name) :: varname !< Variable name.
  integer :: pos !< Domain position.
endtype DomainDimension_t


!> @brief netcdf domain file type.
type, extends(FmsNetcdfFile_t), public :: FmsNetcdfDomainFile_t
  type(domain2d) :: domain !< Two-dimensional domain.
  type(DomainDimension_t), dimension(:), allocatable :: xdims !< Dimensions associated
                                                              !! with the "x" axis
                                                              !! of a 2d domain.
  integer :: nx !< Number of "x" dimensions.
  type(DomainDimension_t), dimension(:), allocatable :: ydims !< Dimensions associated
                                                              !! with the "y" axis
                                                              !! of a 2d domain.
  integer :: ny !< Number of "y" dimensions.
  character(len=256) :: non_mangled_path !< Non-domain-mangled file path.
  logical :: adjust_indices !< Flag telling if indices need to be adjusted
                            !! for domain-decomposed read.
endtype FmsNetcdfDomainFile_t


public :: open_domain_file
public :: close_domain_file
public :: register_domain_decomposed_dimension
public :: register_domain_variable
public :: register_domain_restart_variable_0d
public :: register_domain_restart_variable_1d
public :: register_domain_restart_variable_2d
public :: register_domain_restart_variable_3d
public :: register_domain_restart_variable_4d
public :: register_domain_restart_variable_5d
public :: domain_read_0d
public :: domain_read_1d
public :: domain_read_2d
public :: domain_read_3d
public :: domain_read_4d
public :: domain_read_5d
public :: domain_write_0d
public :: domain_write_1d
public :: domain_write_2d
public :: domain_write_3d
public :: domain_write_4d
public :: domain_write_5d
public :: save_domain_restart
public :: restore_domain_state
public :: get_compute_domain_dimension_indices
public :: get_global_io_domain_indices
public :: is_dimension_registered
public :: get_mosaic_tile_grid

interface compute_global_checksum
  module procedure compute_global_checksum_2d
  module procedure compute_global_checksum_3d
  module procedure compute_global_checksum_4d
end interface compute_global_checksum


contains


!> @brief Get the index of a domain decomposed dimension.
!! @return Index of domain decomposed dimension.
function get_domain_decomposed_index(name_, array, size_) &
  result(index_)

  character(len=*), intent(in) :: name_ !< Name.
  type(DomainDimension_t), dimension(:), intent(in) :: array !< Array to search through.
  integer, intent(in) :: size_ !< Number of spots to look in.
  integer :: index_

  integer :: i

  index_ = variable_not_found
  do i = 1, size_
    if (string_compare(array(i)%varname, name_)) then
      index_ = i
      return
    endif
  enddo
end function get_domain_decomposed_index


!> @brief Add a domain decomposed dimension to an array.
subroutine append_domain_decomposed_dimension(name_, position_, array, size_)

  character(len=*), intent(in) :: name_ !< Variable name.
  integer, intent(in) :: position_ !< Domain position.
  type(DomainDimension_t), dimension(:), intent(inout) :: array !< Array to search through.
  integer, intent(inout) :: size_ !< Number of spots to look in.

  integer :: i

  do i = 1, size_
    if (string_compare(array(i)%varname, name_)) then
      call error("variable "//trim(name_)//" already registered.")
    endif
  enddo
  size_ = size_ + 1
  if (size_ .gt. size(array)) then
    call error("number of domain decomposed variables exceeds limit.")
  endif
  call string_copy(array(size_)%varname, name_)
  array(size_)%pos = position_
end subroutine append_domain_decomposed_dimension


!> @brief Given a domain decomposed dimension, get its domain position.
!! @return Position of the domain decomposed variable.
function get_domain_position(name_, array, size_) &
  result(dpos)

  character(len=*), intent(in) :: name_ !< Variable name.
  type(DomainDimension_t), dimension(:), intent(in) :: array !< Array to search through.
  integer, intent(in) :: size_
  integer :: dpos

  dpos = get_domain_decomposed_index(name_, array, size_)
  if (dpos .ne. variable_not_found) then
    dpos = array(dpos)%pos
  endif
end function get_domain_position


!> @brief Given a variable, get the index of the "x" or "y" domain decomposed
!!        dimension.
!! @return Index of the domain decomposed dimension or else
!!         no_domain_decomposed_dimension.
function get_domain_decomposed_dimension_index(fileobj, variable_name, &
                                               xory, broadcast) &
  result(index_)

  type(FmsNetcdfDomainFile_t), intent(in), target :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  character(len=*), intent(in) :: xory !< String telling which dimension to
                                       !! look for.  Valid values are "x"
                                       !! or "y".
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the index will be
                                             !! broadcast to non "I/O root"
                                             !! ranks.  The broadcast will
                                             !! be done by default.
  integer :: index_

  integer :: ndims
  character(len=nf90_max_name), dimension(:), allocatable :: dim_names
  type(DomainDimension_t), dimension(:), pointer :: p
  integer :: n
  integer :: i

  index_ = no_domain_decomposed_dimension
  if (fileobj%is_root) then
    ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast=.false.)
    allocate(dim_names(ndims))
    dim_names(:) = ""
    call get_variable_dimension_names(fileobj, variable_name, dim_names, broadcast=.false.)
    if (string_compare(xory, x, .true.)) then
      p => fileobj%xdims
      n = fileobj%nx
    elseif (string_compare(xory, y, .true.)) then
      p => fileobj%ydims
      n = fileobj%ny
    else
      call error("unrecognized xory flag value.")
    endif
    do i = 1, size(dim_names)
      if (get_domain_decomposed_index(dim_names(i), p, n) .ne. variable_not_found) then
        index_ = i
        exit
      endif
    enddo
    deallocate(dim_names)
  endif
  if (present(broadcast)) then
    if (.not. broadcast) then
      return
    endif
  endif
  call mpp_broadcast(index_, fileobj%io_root, pelist=fileobj%pelist)
end function get_domain_decomposed_dimension_index


!> @brief Determine if a variable is "domain decomposed."
!! @return Flag telling if the variable is "domain decomposed."
function is_variable_domain_decomposed(fileobj, variable_name, broadcast, &
                                       xindex, yindex, xpos, ypos) &
  result(is_decomposed)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the index will be
                                             !! broadcast to non "I/O root"
                                             !! ranks.  The broadcast will
                                             !! be done by default.
  integer, intent(out), optional :: xindex !< The index of the domain
                                           !! x dimension.
  integer, intent(out), optional :: yindex !< The index of the domain
                                           !! y dimension.
  integer, intent(out), optional :: xpos !< Domain position of the x dimension.
  integer, intent(out), optional :: ypos !< Domain position of the y dimension.
  logical :: is_decomposed

  integer, dimension(2) :: indices
  integer :: ndims
  character(len=nf90_max_name), dimension(:), allocatable :: dim_names

  indices(1) = get_domain_decomposed_dimension_index(fileobj, variable_name, x, broadcast)
  if (present(xindex)) then
    xindex = indices(1)
  endif
  indices(2) = get_domain_decomposed_dimension_index(fileobj, variable_name, y, broadcast)
  if (present(yindex)) then
    yindex = indices(2)
  endif
  is_decomposed = (indices(1) .ne. no_domain_decomposed_dimension) .and. &
                  (indices(2) .ne. no_domain_decomposed_dimension)
  if (is_decomposed) then
    if (.not. present(xpos) .and. .not. present(ypos)) then
      return
    endif
    ndims = get_variable_num_dimensions(fileobj, variable_name, broadcast)
    allocate(dim_names(ndims))
    dim_names(:) = ""
    call get_variable_dimension_names(fileobj, variable_name, dim_names, broadcast)
    if (present(xpos)) then
      xpos = get_domain_position(dim_names(indices(1)), fileobj%xdims, fileobj%nx)
    endif
    if (present(ypos)) then
      ypos = get_domain_position(dim_names(indices(2)), fileobj%ydims, fileobj%ny)
    endif
    deallocate(dim_names)
  else
    if (present(xpos)) then
      xpos = -1
    endif
    if (present(ypos)) then
      ypos = -1
    endif
  endif
end function is_variable_domain_decomposed


!> @brief Determine whether a domain-decomposed dimension has been registered to the file object
!! @return Flag telling if the dimension is registered to the file object
function is_dimension_registered(fileobj, dimension_name) &
  result(is_registered)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: dimension_name !< Dimension name.

  ! local
  logical :: is_registered

  integer :: dpos
  integer :: ndims
  character(len=nf90_max_name), dimension(:), allocatable :: dim_names

  dpos = 0
  is_registered = .false.
  dpos = get_domain_decomposed_index(dimension_name, fileobj%xdims, fileobj%nx)
  if (dpos .ne. variable_not_found) then
    is_registered = .true.
  else
    dpos = get_domain_decomposed_index(dimension_name, fileobj%ydims, fileobj%ny)
    if (dpos .ne. variable_not_found) is_registered = .true.
  endif

end function is_dimension_registered


!> @brief Open a domain netcdf file.
!! @return Flag telling if the open completed successfully.
function open_domain_file(fileobj, path, mode, domain, nc_format, is_restart, dont_add_res_to_filename) &
  result(success)

  type(FmsNetcdfDomainFile_t),intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: path !< File path.
  character(len=*), intent(in) :: mode !< File mode.  Allowed values
                                       !! are "read", "append", "write", or
                                       !! "overwrite".
  type(domain2d), intent(in) :: domain !< Two-dimensional domain.
  character(len=*), intent(in), optional :: nc_format !< Netcdf format that
                                                      !! new files are written
                                                      !! as.  Allowed values
                                                      !! are: "64bit", "classic",
                                                      !! or "netcdf4". Defaults to
                                                      !! "64bit".
  logical, intent(in), optional :: is_restart !< Flag telling if this file
                                              !! is a restart file.  Defaults
                                              !! to false.
  logical, intent(in), optional :: dont_add_res_to_filename !< Flag indicating not to add
                                              !! ".res" to the filename
  logical :: success

  integer, dimension(2) :: io_layout
  integer, dimension(1) :: tile_id
  character(len=256) :: combined_filepath
  type(domain2d), pointer :: io_domain
  character(len=256) :: distributed_filepath
  integer :: pelist_size
  integer, dimension(:), allocatable :: pelist
  logical :: success2
  type(FmsNetcdfDomainFile_t) :: fileobj2

  !Get the path of a "combined" file.
  io_layout = mpp_get_io_domain_layout(domain)
  if (mpp_get_ntile_count(domain) .gt. 1) then
    tile_id = mpp_get_tile_id(domain)
    call domain_tile_filepath_mangle(combined_filepath, path, tile_id(1))
  else
    call string_copy(combined_filepath, path)
  endif

  !Get the path of a "distributed" file.
  io_domain => mpp_get_io_domain(domain)
  if (.not. associated(io_domain)) then
    call error("input domain does not have an io_domain.")
  endif
  if (io_layout(1)*io_layout(2) .gt. 1) then
    tile_id = mpp_get_tile_id(io_domain)
    call io_domain_tile_filepath_mangle(distributed_filepath, combined_filepath, tile_id(1))
  else
    call string_copy(distributed_filepath, combined_filepath)
  endif

  !Make sure the input domain has an I/O domain and get its pelist.
  pelist_size = mpp_get_domain_npes(io_domain)
  allocate(pelist(pelist_size))
  call mpp_get_pelist(io_domain, pelist)
  fileobj%adjust_indices = .true. !Set the default to true

  !Open the distibuted files.
  success = netcdf_file_open(fileobj, distributed_filepath, mode, nc_format, pelist, &
                             is_restart, dont_add_res_to_filename)
  if (string_compare(mode, "read", .true.) .or. string_compare(mode, "append", .true.)) then
    if (success) then
      if (.not. string_compare(distributed_filepath, combined_filepath)) then
        success2 = netcdf_file_open(fileobj2, combined_filepath, mode, nc_format, pelist, &
                                    is_restart, dont_add_res_to_filename)
        if (success2) then
          call error("you have both combined and distributed files.")
        endif
      endif
    else
      success = netcdf_file_open(fileobj, combined_filepath, mode, nc_format, pelist, &
                                 is_restart, dont_add_res_to_filename)
      !If the file is combined and the layout is not (1,1) set the adjust_indices flag to false
      if (success .and. (io_layout(1)*io_layout(2) .gt. 1)) fileobj%adjust_indices = .false.
    endif
  endif
  if (.not. success) then
    deallocate(pelist)
    return
  endif

  !Store/initialize necessary properties.
  call string_copy(fileobj%non_mangled_path, path)
  fileobj%domain = domain
  allocate(fileobj%xdims(max_num_domain_decomposed_dims))
  fileobj%nx = 0
  allocate(fileobj%ydims(max_num_domain_decomposed_dims))
  fileobj%ny = 0
  call string_copy(fileobj%non_mangled_path, path)

  if (string_compare(mode, "write", .true.) .or. string_compare(mode, "overwrite", .true.)) then
    !Add global attribute needed by mppnccombine.
    call register_global_attribute(fileobj, "NumFilesInSet", io_layout(1)*io_layout(2))
  endif
end function open_domain_file


!> @brief Close a domain netcdf file.
subroutine close_domain_file(fileobj)

  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj !< File object.

  call netcdf_file_close(fileobj)
  deallocate(fileobj%xdims)
  fileobj%nx = 0
  deallocate(fileobj%ydims)
  fileobj%ny = 0
end subroutine close_domain_file


!> @brief Add a dimension to a file associated with a two-dimensional domain.
subroutine register_domain_decomposed_dimension(fileobj, dim_name, xory, domain_position)

  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: dim_name !< Dimension name.
  character(len=*), intent(in) :: xory !< Flag telling if the dimension
                                       !! is associated with the "x" or "y"
                                       !! axis of the 2d domain.  Allowed
                                       !! values are "x" or "y".
  integer, intent(in), optional :: domain_position !< Domain position.

  integer :: dpos
  type(domain2d), pointer :: io_domain
  integer :: domain_size
  integer :: dim_size

  dpos = default_domain_position
  if (mpp_domain_is_symmetry(fileobj%domain) .and. present(domain_position)) then
    dpos = domain_position
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  if (string_compare(xory, x, .true.)) then
    if (dpos .ne. center .and. dpos .ne. east) then
      call error("only center or east supported for x dimensions.")
    endif
    call mpp_get_global_domain(io_domain, xsize=domain_size, position=dpos)
    call append_domain_decomposed_dimension(dim_name, dpos, fileobj%xdims, fileobj%nx)
  elseif (string_compare(xory, y, .true.)) then
    if (dpos .ne. center .and. dpos .ne. north) then
      call error("only center or north supported for y dimensions.")
    endif
    call mpp_get_global_domain(io_domain, ysize=domain_size, position=dpos)
    call append_domain_decomposed_dimension(dim_name, dpos, fileobj%ydims, fileobj%ny)
  else
    call error("unrecognized xory flag value.")
  endif
  if (fileobj%is_readonly .or. (fileobj%mode_is_append .and. dimension_exists(fileobj, dim_name))) then
    call get_dimension_size(fileobj, dim_name, dim_size, broadcast=.true.)
    if (dim_size .lt. domain_size) then
      call error("dimension "//trim(dim_name)//" is smaller than the size of" &
                 //" the associated domain "//trim(xory)//" axis.")
    endif
  else
    call netcdf_add_dimension(fileobj, dim_name, domain_size, is_compressed=.false.)
  endif
end subroutine register_domain_decomposed_dimension


!> @brief Add a "domain_decomposed" attribute to certain variables because it is
!!        required by mppnccombine.
!! @internal
subroutine add_domain_attribute(fileobj, variable_name)

  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.

  type(domain2d), pointer :: io_domain
  integer :: dpos
  integer :: sg
  integer :: eg
  integer :: s
  integer :: e

  io_domain => mpp_get_io_domain(fileobj%domain)
  dpos = get_domain_decomposed_index(variable_name, fileobj%xdims, fileobj%nx)
  if (dpos .ne. variable_not_found) then
    dpos = fileobj%xdims(dpos)%pos
    call mpp_get_global_domain(fileobj%domain, xbegin=sg, xend=eg, position=dpos)
    call mpp_get_global_domain(io_domain, xbegin=s, xend=e, position=dpos)
    call register_variable_attribute(fileobj, variable_name, "domain_decomposition", &
                                     (/sg, eg, s, e/))
  else
    dpos = get_domain_decomposed_index(variable_name, fileobj%ydims, fileobj%ny)
    if (dpos .ne. variable_not_found) then
      dpos = fileobj%ydims(dpos)%pos
      call mpp_get_global_domain(fileobj%domain, ybegin=sg, yend=eg, position=dpos)
      call mpp_get_global_domain(io_domain, ybegin=s, yend=e, position=dpos)
      call register_variable_attribute(fileobj, variable_name, "domain_decomposition", &
                                       (/sg, eg, s, e/))
    endif
  endif
end subroutine add_domain_attribute


!> @brief Add a domain decomposed variable.
subroutine register_domain_variable(fileobj, variable_name, variable_type, dimensions)

  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  character(len=*), intent(in) :: variable_type !< Variable type.  Allowed
                                                !! values are: "int", "int64",
                                                !! "float", or "double".
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.

  if (.not. fileobj%is_readonly) then
    call netcdf_add_variable(fileobj, variable_name, variable_type, dimensions)
    if (present(dimensions)) then
      if (size(dimensions) .eq. 1) then
        call add_domain_attribute(fileobj, variable_name)
      endif
    endif
  endif
end subroutine register_domain_variable


!> @brief Loop through registered restart variables and write them to
!!        a netcdf file.
subroutine save_domain_restart(fileobj, unlim_dim_level)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.

  integer :: i
  character(len=32) :: chksum
  logical :: is_decomposed

  if (.not. fileobj%is_restart) then
    call error("file "//trim(fileobj%path)//" is not a restart file.")
  endif

! Calculate the variable's checksum and write it to the netcdf file
  do i = 1, fileobj%num_restart_vars
    if (associated(fileobj%restart_vars(i)%data2d)) then
      chksum = compute_global_checksum(fileobj, fileobj%restart_vars(i)%varname, &
                                       fileobj%restart_vars(i)%data2d, is_decomposed)
      if (is_decomposed) then
        call register_variable_attribute(fileobj, fileobj%restart_vars(i)%varname, &
                                         "checksum", chksum, str_len=len(chksum))
      endif
    elseif (associated(fileobj%restart_vars(i)%data3d)) then
      chksum = compute_global_checksum(fileobj, fileobj%restart_vars(i)%varname, &
                                       fileobj%restart_vars(i)%data3d, is_decomposed)
      if (is_decomposed) then
        call register_variable_attribute(fileobj, fileobj%restart_vars(i)%varname, &
                                         "checksum", chksum, str_len=len(chksum))
      endif
    elseif (associated(fileobj%restart_vars(i)%data4d)) then
      chksum = compute_global_checksum(fileobj, fileobj%restart_vars(i)%varname, &
                                       fileobj%restart_vars(i)%data4d, is_decomposed)
      if (is_decomposed) then
        call register_variable_attribute(fileobj, fileobj%restart_vars(i)%varname, &
                                         "checksum", chksum, str_len=len(chksum))
      endif
    endif
  enddo

! Write the variable's data to the netcdf file
  do i = 1, fileobj%num_restart_vars
    if (associated(fileobj%restart_vars(i)%data0d)) then
      call domain_write_0d(fileobj, fileobj%restart_vars(i)%varname, &
                           fileobj%restart_vars(i)%data0d, unlim_dim_level=unlim_dim_level)
    elseif (associated(fileobj%restart_vars(i)%data1d)) then
      call domain_write_1d(fileobj, fileobj%restart_vars(i)%varname, &
                           fileobj%restart_vars(i)%data1d, unlim_dim_level=unlim_dim_level)
    elseif (associated(fileobj%restart_vars(i)%data2d)) then
      call domain_write_2d(fileobj, fileobj%restart_vars(i)%varname, &
                           fileobj%restart_vars(i)%data2d, unlim_dim_level=unlim_dim_level)
    elseif (associated(fileobj%restart_vars(i)%data3d)) then
      call domain_write_3d(fileobj, fileobj%restart_vars(i)%varname, &
                           fileobj%restart_vars(i)%data3d, unlim_dim_level=unlim_dim_level)
    elseif (associated(fileobj%restart_vars(i)%data4d)) then
      call domain_write_4d(fileobj, fileobj%restart_vars(i)%varname, &
                           fileobj%restart_vars(i)%data4d, unlim_dim_level=unlim_dim_level)
    else
      call error("This routine only accepts data that is scalar, 1d 2d 3d or 4d.  The data sent in has an unsupported dimensionality")
    endif
  enddo

end subroutine save_domain_restart


!> @brief Loop through registered restart variables and read them from
!!        a netcdf file.
subroutine restore_domain_state(fileobj, unlim_dim_level)

  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj !< File object.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.

  integer :: i
  character(len=32) :: chksum_in_file
  character(len=32) :: chksum
  logical :: is_decomposed

  if (.not. fileobj%is_restart) then
    call error("file "//trim(fileobj%path)//" is not a restart file.")
  endif
  do i = 1, fileobj%num_restart_vars
    if (associated(fileobj%restart_vars(i)%data0d)) then
      call domain_read_0d(fileobj, fileobj%restart_vars(i)%varname, &
                          fileobj%restart_vars(i)%data0d, unlim_dim_level=unlim_dim_level)
    elseif (associated(fileobj%restart_vars(i)%data1d)) then
      call domain_read_1d(fileobj, fileobj%restart_vars(i)%varname, &
                          fileobj%restart_vars(i)%data1d, unlim_dim_level=unlim_dim_level)
    elseif (associated(fileobj%restart_vars(i)%data2d)) then
      call domain_read_2d(fileobj, fileobj%restart_vars(i)%varname, &
                          fileobj%restart_vars(i)%data2d, unlim_dim_level=unlim_dim_level)
      chksum = compute_global_checksum(fileobj, fileobj%restart_vars(i)%varname, &
                                       fileobj%restart_vars(i)%data2d, is_decomposed)
      if (variable_att_exists(fileobj, fileobj%restart_vars(i)%varname, "checksum") .and. &
          is_decomposed) then
        call get_variable_attribute(fileobj, fileobj%restart_vars(i)%varname, &
                                    "checksum", chksum_in_file)
        if (.not. string_compare(trim(adjustl(chksum_in_file)), trim(adjustl(chksum)))) then
          call error("checksum attribute does not match data in file.")
        endif
      endif
    elseif (associated(fileobj%restart_vars(i)%data3d)) then
      call domain_read_3d(fileobj, fileobj%restart_vars(i)%varname, &
                          fileobj%restart_vars(i)%data3d, unlim_dim_level=unlim_dim_level)
      chksum = compute_global_checksum(fileobj, fileobj%restart_vars(i)%varname, &
                                       fileobj%restart_vars(i)%data3d, is_decomposed)
      if (variable_att_exists(fileobj, fileobj%restart_vars(i)%varname, "checksum") .and. &
          is_decomposed) then
        call get_variable_attribute(fileobj, fileobj%restart_vars(i)%varname, &
                                    "checksum", chksum_in_file)
        if (.not. string_compare(trim(adjustl(chksum_in_file)), trim(adjustl(chksum)))) then
          call error("checksum attribute does not match data in file.")
        endif
      endif
    elseif (associated(fileobj%restart_vars(i)%data4d)) then
      call domain_read_4d(fileobj, fileobj%restart_vars(i)%varname, &
                          fileobj%restart_vars(i)%data4d, unlim_dim_level=unlim_dim_level)
      chksum = compute_global_checksum(fileobj, fileobj%restart_vars(i)%varname, &
                                       fileobj%restart_vars(i)%data4d, is_decomposed)
      if (variable_att_exists(fileobj, fileobj%restart_vars(i)%varname, "checksum") .and. &
          is_decomposed) then
        call get_variable_attribute(fileobj, fileobj%restart_vars(i)%varname, &
                                    "checksum", chksum_in_file)
        if (.not. string_compare(trim(adjustl(chksum_in_file)), trim(adjustl(chksum)))) then
          call error("checksum attribute does not match data in file.")
        endif
      endif
    else
      call error("this branch should not be reached.")
    endif
  enddo
end subroutine restore_domain_state


!> @brief Return an array of compute domain indices.
subroutine get_compute_domain_dimension_indices(fileobj, dimname, indices)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: dimname !< Name of dimension variable.
  integer, dimension(:), allocatable, intent(inout) :: indices !< Compute domain indices.

  type(domain2d), pointer :: io_domain
  integer :: dpos
  integer :: s
  integer :: e
  integer :: i

  io_domain => mpp_get_io_domain(fileobj%domain)
  dpos = get_domain_decomposed_index(dimname, fileobj%xdims, fileobj%nx)
  if (dpos .ne. variable_not_found) then
    dpos = fileobj%xdims(dpos)%pos
    call mpp_get_compute_domain(io_domain, xbegin=s, xend=e, position=dpos)
  else
    dpos = get_domain_decomposed_index(dimname, fileobj%ydims, fileobj%ny)
    if (dpos .ne. variable_not_found) then
      dpos = fileobj%ydims(dpos)%pos
      call mpp_get_compute_domain(io_domain, ybegin=s, yend=e, position=dpos)
    else
      call error("input dimension is not associated with the domain.")
    endif
  endif
  if (allocated(indices)) then
    deallocate(indices)
  endif
  allocate(indices(e-s+1))
  do i = s, e
    indices(i-s+1) = i
  enddo
end subroutine get_compute_domain_dimension_indices


!> @brief Utility routine that retrieves domain indices.
!! @internal
subroutine domain_offsets(data_xsize, data_ysize, domain, xpos, ypos, &
                          isd, isc, xc_size, jsd, jsc, yc_size, &
                          buffer_includes_halos, extra_x_point, &
                          extra_y_point)

  integer, intent(in) :: data_xsize !< Size of buffer's domain "x" dimension.
  integer, intent(in) :: data_ysize !< Size of buffer's domain "y" dimension.
  type(domain2d), intent(in) :: domain !< Parent domain.
  integer, intent(in) :: xpos !< Variable's domain x dimension position.
  integer, intent(in) :: ypos !< Variable's domain y dimension position.
  integer, intent(out) :: isd !< Starting index for x dimension of data domain.
  integer, intent(out) :: isc !< Starting index for x dimension of compute domain.
  integer, intent(out) :: xc_size !< Size of x dimension of compute domain.
  integer, intent(out) :: jsd !< Starting index for y dimension of data domain.
  integer, intent(out) :: jsc !< Starting index for y dimension of compute domain.
  integer, intent(out) :: yc_size !< Size of y dimension of compute domain.
  logical, intent(out) :: buffer_includes_halos !< Flag telling if input buffer includes space for halos.
  logical, intent(out), optional :: extra_x_point !<
  logical, intent(out), optional :: extra_y_point !<

  integer :: xd_size
  integer :: yd_size
  integer :: iec
  integer :: xmax
  integer :: jec
  integer :: ymax
  type(domain2d), pointer :: io_domain !< I/O domain variable is decomposed over.

  io_domain => mpp_get_io_domain(domain)

  call mpp_get_global_domain(domain, xend=xmax, position=xpos)
  call mpp_get_global_domain(domain, yend=ymax, position=ypos)
  call mpp_get_data_domain(io_domain, xbegin=isd, xsize=xd_size, position=xpos)
  call mpp_get_data_domain(io_domain, ybegin=jsd, ysize=yd_size, position=ypos)

  call mpp_get_compute_domain(io_domain, xbegin=isc, xend=iec, xsize=xc_size, &
                              position=xpos)
  ! If the xpos is east and the ending x index is NOT equal to max allowed, set extra_x_point to true
  if (present(extra_x_point)) then
    if ((xpos .eq. east) .and. (iec .ne. xmax)) then
      extra_x_point = .true.
    else
      extra_x_point = .false.
    endif
  endif

  call mpp_get_compute_domain(io_domain, ybegin=jsc, yend=jec, ysize=yc_size, &
                              position=ypos)
  ! If the ypost is north and the ending y index is NOT equal to max allowed, set extra_y_point to true
  if (present(extra_y_point)) then
    if ((ypos .eq. north) .and. (jec .ne. ymax)) then
      extra_y_point = .true.
    else
      extra_y_point = .false.
    endif
  endif

  buffer_includes_halos = (data_xsize .eq. xd_size) .and. (data_ysize .eq. yd_size)
  if (.not. buffer_includes_halos .and. data_xsize .ne. xc_size .and. data_ysize &
      .ne. yc_size) then
    call error("size of x dimension of input buffer does not match size" &
               //" of x dimension of data or compute domain.")
  endif
end subroutine domain_offsets


!> @brief Get starting/ending global indices of the I/O domain for a domain decomposed
!!        file.
subroutine get_global_io_domain_indices(fileobj, dimname, is, ie, indices)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: dimname !< Name of dimension variable.
  integer, intent(out) :: is !< Staring index of I/O global domain.
  integer, intent(out) :: ie !< Ending index of I/O global domain.
  integer, dimension(:), allocatable, intent(out), optional :: indices !< Global domain indices

  type(domain2d), pointer :: io_domain
  integer :: dpos
  integer :: i

  io_domain => mpp_get_io_domain(fileobj%domain)
  dpos = get_domain_decomposed_index(dimname, fileobj%xdims, fileobj%nx)
  if (dpos .ne. variable_not_found) then
    dpos = fileobj%xdims(dpos)%pos
    call mpp_get_global_domain(io_domain, xbegin=is, xend=ie, position=dpos)
  else
    dpos = get_domain_decomposed_index(dimname, fileobj%ydims, fileobj%ny)
    if (dpos .ne. variable_not_found) then
      dpos = fileobj%ydims(dpos)%pos
      call mpp_get_global_domain(io_domain, ybegin=is, yend=ie, position=dpos)
    else
      call error("input dimension is not associated with the domain.")
    endif
  endif

! Allocate indices to the difference between the ending and starting indices and
! fill indices with the data
  if (present(indices)) then
    if(allocated(indices)) then
      call error("get_global_io_domain_indices: the variable indices should not be allocated.")
    endif
    allocate(indices(ie-is+1))
    do i = is, ie
      indices(i-is+1) = i
    enddo
  endif


end subroutine get_global_io_domain_indices

!> @brief Read a mosaic_file and get the grid filename for the current tile or
!!        for the tile specified
subroutine get_mosaic_tile_grid(grid_file,mosaic_file, domain, tile_count)
  character(len=*), intent(out)          :: grid_file !< Filename of the grid file for the
                                                      !! current domain tile or for tile
                                                      !! specified in tile_count
  character(len=*), intent(in)           :: mosaic_file !< Filename that will be read
  type(domain2D),   intent(in)           :: domain !< Input domain
  integer,          intent(in), optional :: tile_count !< Optional argument indicating
                                                       !! the tile you want grid file name for
                                                       !! this is for when a pe is in more than
                                                       !! tile.
  integer                                :: tile !< Current domian tile or tile_count
  integer                                :: ntileMe !< Total number of tiles in the domain
  integer, dimension(:), allocatable     :: tile_id !< List of tiles in the domain
  type(FmsNetcdfFile_t)                  :: fileobj !< Fms2io file object

  tile = 1
  if(present(tile_count)) tile = tile_count
  ntileMe = mpp_get_current_ntile(domain)
  allocate(tile_id(ntileMe))
  tile_id = mpp_get_tile_id(domain)

  if (netcdf_file_open(fileobj, mosaic_file, "read")) then
      call netcdf_read_data(fileobj, "gridfiles", grid_file, corner=tile_id(tile))
      grid_file = 'INPUT/'//trim(grid_file)
      call netcdf_file_close(fileobj)
  endif

end subroutine get_mosaic_tile_grid

include "register_domain_restart_variable.inc"
include "domain_read.inc"
include "domain_write.inc"
#include "compute_global_checksum.inc"


end module fms_netcdf_domain_io_mod

