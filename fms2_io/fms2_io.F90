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

!> @brief Public API for fms I/O.
module fms2_io_mod
use fms_io_utils_mod
use netcdf_io_mod
use fms_netcdf_domain_io_mod
use fms_netcdf_unstructured_domain_io_mod
use blackboxio
use mpp_mod, only: mpp_init, input_nml_file, mpp_error, FATAL
use mpp_domains_mod, only: mpp_domains_init
implicit none
private


public :: unlimited
public :: FmsNetcdfFile_t
public :: FmsNetcdfDomainFile_t
public :: FmsNetcdfUnstructuredDomainFile_t
public :: open_file
public :: open_virtual_file
public :: close_file
public :: register_axis
public :: register_field
public :: register_restart_field
public :: write_data
public :: read_data
public :: write_restart
public :: write_new_restart
public :: read_restart
public :: read_new_restart
public :: global_att_exists
public :: variable_att_exists
public :: register_global_attribute
public :: register_variable_attribute
public :: get_global_attribute
public :: get_variable_attribute
public :: get_num_dimensions
public :: get_dimension_names
public :: dimension_exists
public :: is_dimension_unlimited
public :: get_dimension_size
public :: get_num_variables
public :: get_variable_names
public :: variable_exists
public :: get_variable_num_dimensions
public :: get_variable_dimension_names
public :: get_variable_size
public :: get_compute_domain_dimension_indices
public :: get_global_io_domain_indices
public :: Valid_t
public :: get_valid
public :: is_valid
public :: get_unlimited_dimension_name
public :: get_variable_unlimited_dimension_index
public :: file_exists
public :: compressed_start_and_count
public :: get_variable_sense
public :: get_variable_missing
public :: get_variable_units
public :: get_time_calendar
public :: open_check
public :: is_registered_to_restart
public :: check_if_open
public :: set_fileobj_time_name
public :: is_dimension_registered
public :: fms2_io_init
public :: get_mosaic_tile_grid

interface open_file
  module procedure netcdf_file_open_wrap
  module procedure open_domain_file
  module procedure open_unstructured_domain_file
end interface open_file


interface open_virtual_file
  module procedure create_diskless_netcdf_file_wrap
  module procedure create_diskless_domain_file
  module procedure create_diskless_unstructured_domain_file
end interface open_virtual_file


interface close_file
  module procedure netcdf_file_close_wrap
  module procedure close_domain_file
  module procedure close_unstructured_domain_file
end interface close_file


interface register_axis
  module procedure netcdf_add_dimension
  module procedure register_compressed_dimension
  module procedure register_domain_decomposed_dimension
  module procedure register_unstructured_dimension
end interface register_axis


interface register_field
  module procedure netcdf_add_variable_wrap
  module procedure register_domain_variable
  module procedure register_unstructured_domain_variable
end interface register_field


interface register_restart_field
  module procedure netcdf_add_restart_variable_0d_wrap
  module procedure netcdf_add_restart_variable_1d_wrap
  module procedure netcdf_add_restart_variable_2d_wrap
  module procedure netcdf_add_restart_variable_3d_wrap
  module procedure netcdf_add_restart_variable_4d_wrap
  module procedure netcdf_add_restart_variable_5d_wrap
  module procedure register_domain_restart_variable_0d
  module procedure register_domain_restart_variable_1d
  module procedure register_domain_restart_variable_2d
  module procedure register_domain_restart_variable_3d
  module procedure register_domain_restart_variable_4d
  module procedure register_domain_restart_variable_5d
  module procedure register_unstructured_domain_restart_variable_0d
  module procedure register_unstructured_domain_restart_variable_1d
  module procedure register_unstructured_domain_restart_variable_2d
  module procedure register_unstructured_domain_restart_variable_3d
  module procedure register_unstructured_domain_restart_variable_4d
  module procedure register_unstructured_domain_restart_variable_5d
end interface register_restart_field


interface write_data
  module procedure compressed_write_0d_wrap
  module procedure compressed_write_1d_wrap
  module procedure compressed_write_2d_wrap
  module procedure compressed_write_3d_wrap
  module procedure compressed_write_4d_wrap
  module procedure compressed_write_5d_wrap
  module procedure domain_write_0d
  module procedure domain_write_1d
  module procedure domain_write_2d
  module procedure domain_write_3d
  module procedure domain_write_4d
  module procedure domain_write_5d
  module procedure unstructured_domain_write_0d
  module procedure unstructured_domain_write_1d
  module procedure unstructured_domain_write_2d
  module procedure unstructured_domain_write_3d
  module procedure unstructured_domain_write_4d
  module procedure unstructured_domain_write_5d
end interface write_data


interface read_data
  module procedure compressed_read_0d
  module procedure compressed_read_1d
  module procedure compressed_read_2d
  module procedure compressed_read_3d
  module procedure compressed_read_4d
  module procedure compressed_read_5d
  module procedure domain_read_0d
  module procedure domain_read_1d
  module procedure domain_read_2d
  module procedure domain_read_3d
  module procedure domain_read_4d
  module procedure domain_read_5d
  module procedure unstructured_domain_read_0d
  module procedure unstructured_domain_read_1d
  module procedure unstructured_domain_read_2d
  module procedure unstructured_domain_read_3d
  module procedure unstructured_domain_read_4d
  module procedure unstructured_domain_read_5d
end interface read_data


interface write_restart
  module procedure netcdf_save_restart_wrap
  module procedure save_domain_restart
  module procedure unstructured_write_restart
end interface write_restart


interface write_new_restart
  module procedure netcdf_save_restart_wrap2
  module procedure save_domain_restart_wrap
  module procedure unstructured_write_restart_wrap
end interface write_new_restart


interface read_restart
  module procedure netcdf_restore_state
  module procedure restore_domain_state
end interface read_restart


interface read_new_restart
  module procedure netcdf_restore_state_wrap
  module procedure restore_domain_state_wrap
end interface read_new_restart

logical, private :: fms2_io_is_initialized = .false. !< True after fms2_io_init is run
!< Namelist variables
integer :: ncchksz = 64*1024  !< User defined chunksize (in bytes) argument in netcdf file
                              !! creation calls. Replaces setting the NC_CHKSZ environment variable.
character (len = 10) :: netcdf_default_format = "64bit" !< User defined netcdf file format, acceptable values
                              !! are: "64bit", "classic", "netcdf4". This can be overwritten if you specify
                              !! "nc_format" in the open_file call
integer :: header_buffer_val = 16384 !< Use defined netCDF header buffer size(in bytes) used in
                                     !! NF__ENDDEF
namelist / fms2_io_nml / &
                      ncchksz, netcdf_default_format, header_buffer_val

contains

!> @brief Reads the fms2_io_nml
subroutine fms2_io_init ()
 integer :: mystat

!> Check if the module has already been initialized
  if (fms2_io_is_initialized) return
!> Call initialization routines that this module depends on
  call mpp_init()
  call mpp_domains_init()
!> Read the namelist
  READ (input_nml_file, NML=fms2_io_nml, IOSTAT=mystat)
!>Send the namelist variables to their respective modules
  if (ncchksz .le. 0) then
        call mpp_error(FATAL, "ncchksz in fms2_io_nml must be a positive number.")
  endif
  call netcdf_io_init (ncchksz,header_buffer_val,netcdf_default_format)
  if (header_buffer_val .le. 0) then
        call mpp_error(FATAL, "header_buffer_val in fms2_io_nml must be a positive number.")
  endif
  call blackboxio_init (ncchksz)
!> Mark the fms2_io as initialized
  fms2_io_is_initialized = .true.
end subroutine fms2_io_init


end module fms2_io_mod
