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
module fms_netcdf_unstructured_domain_io_mod
use,intrinsic :: iso_fortran_env
use netcdf
use mpp_domains_mod
use fms_io_utils_mod
use netcdf_io_mod
implicit none
private


!> @brief netcdf unstructured domain file type.
type, extends(FmsNetcdfFile_t), public :: FmsNetcdfUnstructuredDomainFile_t
  type(domainug) :: domain !< Unstructured domain.
  character(len=256) :: non_mangled_path !< Non-domain-mangled path.
endtype FmsNetcdfUnstructuredDomainFile_t


public :: open_unstructured_domain_file
public :: close_unstructured_domain_file
public :: register_unstructured_dimension
public :: register_unstructured_domain_variable
public :: register_unstructured_domain_restart_variable_0d
public :: register_unstructured_domain_restart_variable_1d
public :: register_unstructured_domain_restart_variable_2d
public :: register_unstructured_domain_restart_variable_3d
public :: register_unstructured_domain_restart_variable_4d
public :: register_unstructured_domain_restart_variable_5d
public :: unstructured_domain_read_0d
public :: unstructured_domain_read_1d
public :: unstructured_domain_read_2d
public :: unstructured_domain_read_3d
public :: unstructured_domain_read_4d
public :: unstructured_domain_read_5d
public :: unstructured_domain_write_0d
public :: unstructured_domain_write_1d
public :: unstructured_domain_write_2d
public :: unstructured_domain_write_3d
public :: unstructured_domain_write_4d
public :: unstructured_domain_write_5d
public :: unstructured_write_restart


contains


!> @brief Open a netcdf file that is associated with an unstructured domain.
!! @return Flag telling if the open completed successfully.
function open_unstructured_domain_file(fileobj, path, mode, domain, nc_format, &
                                       is_restart, dont_add_res_to_filename) &
  result(success)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: path !< File path.
  character(len=*), intent(in) :: mode !< File mode.  Allowed values
                                       !! are "read", "append", "write", or
                                       !! "overwrite".
  type(domainug), intent(in) :: domain !< Unstructured domain.
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

  type(domainug), pointer :: io_domain
  integer :: pelist_size
  integer, dimension(:), allocatable :: pelist
  character(len=256) :: buf
  character(len=256) :: buf2
  integer :: tile_id

  !Get the input domain's I/O domain pelist.
  io_domain => mpp_get_ug_io_domain(domain)
  if (.not. associated(io_domain)) then
    call error("input domain does not have an io_domain.")
  endif
  pelist_size = mpp_get_ug_domain_npes(io_domain)
  allocate(pelist(pelist_size))
  call mpp_get_ug_domain_pelist(io_domain, pelist)

  !Add the domain tile id to the file name (if necessary).
  call string_copy(buf, path)
  if (mpp_get_UG_domain_ntiles(domain) .gt. 1) then
    tile_id = mpp_get_ug_domain_tile_id(domain)
    call domain_tile_filepath_mangle(buf, path, tile_id)
  endif

  success = .false.
  if (string_compare(mode, "read", .true.) .or. string_compare(mode, "append", .true.)) then
    !Only for reading: attempt to open non-distributed files.
    success = netcdf_file_open(fileobj, buf, mode, nc_format, pelist, is_restart, dont_add_res_to_filename)
  endif
  if (.not. success) then
    !Add the domain tile id to the file name (if necessary).
    if (mpp_get_io_domain_ug_layout(domain) .gt. 1) then
      tile_id = mpp_get_ug_domain_tile_id(io_domain)
      call string_copy(buf2, buf)
      call io_domain_tile_filepath_mangle(buf, buf2, tile_id)
    endif

    !Open distributed files.
    success = netcdf_file_open(fileobj, buf, mode, nc_format, pelist, is_restart, dont_add_res_to_filename)
  endif
  deallocate(pelist)

  if (.not. success) then
    !This branch should only be entered if the file attempting to be read
    !does not exist.
    return
  endif

  !Store/initialize necessary properties.
  fileobj%domain = domain
  call string_copy(fileobj%non_mangled_path, path)
end function open_unstructured_domain_file


!> @brief Wrapper to distinguish interfaces.
subroutine close_unstructured_domain_file(fileobj)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(inout) :: fileobj !< File object.

  call netcdf_file_close(fileobj)
end subroutine close_unstructured_domain_file


!> @brief Add an unstructured dimension.
subroutine register_unstructured_dimension(fileobj, dim_name)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: dim_name !< Dimension name.

  type(domainug),pointer :: io_domain
  integer, dimension(:), allocatable :: c
  integer, dimension(:), allocatable :: e

  allocate(c(size(fileobj%pelist)))
  allocate(e(size(fileobj%pelist)))
  io_domain => mpp_get_ug_io_domain(fileobj%domain)
  call mpp_get_ug_compute_domains(io_domain, begin=c, size=e)
  if (c(1) .ne. 1) then
    c(:) = c(:) - c(1) + 1
  endif
  call register_compressed_dimension(fileobj, dim_name, c, e)
  deallocate(c)
  deallocate(e)
end subroutine register_unstructured_dimension


!> @brief Wrapper to distinguish interfaces.
subroutine register_unstructured_domain_variable(fileobj, variable_name, &
                                                 variable_type, dimensions)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  character(len=*), intent(in) :: variable_type !< Variable type.  Allowed
                                                !! values are: "int", "int64",
                                                !! "float", or "double".
  character(len=*), dimension(:), intent(in), optional :: dimensions !< Dimension names.

  call netcdf_add_variable(fileobj, variable_name, variable_type, dimensions)
end subroutine register_unstructured_domain_variable


!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_write_restart(fileobj, unlim_dim_level)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.

  call netcdf_save_restart(fileobj, unlim_dim_level)
end subroutine unstructured_write_restart


include "register_unstructured_domain_restart_variable.inc"
include "unstructured_domain_read.inc"
include "unstructured_domain_write.inc"


end module fms_netcdf_unstructured_domain_io_mod
