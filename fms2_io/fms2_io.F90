!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> @defgroup fms2_io_mod fms2_io_mod
!> @ingroup fms2_io
!> @brief This module aims support netCDF I/O operations for 3 use-cases within the FMS framework:
!!
!! 1) Netcdf operations done on a single core via the netcdf_io_mod module.
!! 2) Netcdf operations done with multiple cores on a structured grid via the fms_netcdf_domain_io_mod module.
!! 3) Netcdf operations done with multiple cores on a unstructured grid via the fms_netcdf_unstructured_domain_io_mod
!!    module.
!!
!! fms2_io_mod is the top level module that provides open, close, read, and write interfaces to the NetCDF package.
!! This module defines public "aliases"(interfaces) to select procedures in fms_netcdf_domain_io_mod for reading/writing
!! structured grid domains; fms_netcdf_unstructured_domain_io_mod for reading/writing data on unstructured grid domains;
!! netcdf_io_mod for reading/writing netcdf files on a single core. Subroutines and functions in the latter
!! mentioned modules (fms_netcdf_domain_io_mod, fms_netcdf_unstructured_domain_io_mod, and netcdf_io_mod) are intended
!! for internally use only. We highly recommended to only call public interfaces defined in this module.
!!
!! Before any fms2_io_mod I/O operations, a file derived type must be declared first.
!! Three file derived types are currently available and are described below.  Any instances
!! of these three file types are referred to as "fileobj" in this module.
!!
!! - FmsNetcdfFile_t: provides limited number of wrapper procedures to the netCDF4 library.  If the
!! user provides a pelist to procedures compatible with this type, only the root rank of the pelist
!! performs I/O operations by calling the NetCDF library: the root rank either boradcasts the read-in
!! data to the remaining ranks in the pelist, or gathers data from the remaining ranks in the pelist before writing.
!!
!! - FmsNetcdfDomainFile_t:  extends upon FmsNetcdfFile_t and adds supports for "domain-decomposed" reads and writes.
!!   Here, "domain decomposed" refers to data that is partitioned into subdomains of the decomposed global domain,
!!   and each MPI rank holds its portion of the global data.  The users must provide a domain of type Domain2D from
!!   mpp_domains_mod when initializing this file object.
!!
!! - FmsNetcdfUnstructuredDomainFile_t: extends upon FmsNetcdfUnstructuredDomainFile_t and adds support for
!!   “domain-decomposed” reads/writes for data decomposed on subdomains of unstructured grids
!!   The users must provide a domain of type DomainUG from mpp_domains_mod when initializing this file object.
!!
!! See mpp_domains_mod documentation for more information on creating a domain decomposition using FMS.
!!
!! When doing IO on multiple cores, either via the FmsNetcdfDomainFile_t or the FmsNetcdfUnstructuredDomainFile_t,
!! the number of cores performing IO operations is set by the io_layout. The io_layout is 2 integers set via
!! mpp_set_io_domain and determines how to divide up the existing domain decomposition layout for IO operations.
!! See FmsNetcdfDomainFile_t for more information on how the io_layout works.
!!
!! Besides standard open/close/read/write operations, this module also provides interfaces for writing and reading
!! "diskless" netcdf files via the blackboxio module.
!!
!! Users can specify I/O specifications with the fms2_io_nml namelist which allows users to specify the following
!! NetCDF library parameters:  Netcdf file format, chunk_size, deflate_level, shuffle.  Note, fms2_io accepts the
!! Netcdf file format namelist values:"64bit", "class", and "netcdf4". For more information on optimizing NetCDF writing
!! operation, see https://docs.unidata.ucar.edu/netcdf-c/current/file_format_specifications.html
!!
!! @note The legacy IO modules, fms_io_mod and mpp_io_mod, are
!! If converting legacy code from fms_io/mpp_io to fms2_io, please refer to the migration guide at fms2_io/readme.md.

!> @addtogroup fms2_io_mod
!> @{
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

!> NetCDF constant (enum) for unlimited dimension identification
public :: unlimited

!> File object types are defined in the helper modules (netcdf_io_mod,fms_netcdf_domain_io_mod,
!! fms_netcdf_unstructured_domain_io_mod) but are made public here for user access.
public :: FmsNetcdfFile_t, FmsNetcdfDomainFile_t, FmsNetcdfUnstructuredDomainFile_t

!> Interfaces defined below to make public
public :: open_file, open_virtual_file, close_file
public :: register_axis
public :: register_field
public :: register_restart_field
public :: write_data
public :: read_data
public :: write_restart
public :: write_new_restart
public :: read_restart
public :: read_new_restart

!> Routines/functions from netcdf_io_mod to make public
public :: register_unlimited_compressed_axis
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
public :: flush_file
public :: write_restart_bc
public :: read_restart_bc
public :: compressed_start_and_count
public :: get_variable_sense
public :: get_variable_missing
public :: get_variable_units
public :: get_time_calendar
public :: is_registered_to_restart
public :: check_if_open
public :: set_fileobj_time_name
public :: Valid_t
public :: get_valid
public :: is_valid

!> Routines/functions from fms_netcdf_domain_io_mod to make public
public :: get_compute_domain_dimension_indices
public :: get_global_io_domain_indices
public :: get_unlimited_dimension_name
public :: get_variable_unlimited_dimension_index

!> Routines/functions from fms_io_utils_mod to make public
public :: file_exists
public :: open_check
public :: is_dimension_registered
public :: fms2_io_init
public :: get_mosaic_tile_grid
public :: ascii_read
public :: get_mosaic_tile_file
public :: parse_mask_table
public :: get_filename_appendix
public :: set_filename_appendix
public :: get_instance_filename
public :: nullify_filename_appendix
!> @}

!> @brief Opens a netcdf dataset on disk. Initializes the file object.
!!
!> <br>Example usage:
!!
!!              io_success = open_file(fileobj, "filename", "write")
!!
!! Opens a netcdf file for a standard data, domain decomposed data, or unstructured domain decomposed data based off
!! the fileobj used, and also intializes the fileobj for subsequent IO operations.
!!
!! File mode is set to one of "read"/"write"/"overwrite"/"append":
!!
!!              io_success = open_file(fileobj, "filename", "overwrite", domain)
!!
!! Opens a domain netcdf file of type @ref fmsnetcdfdomainfile_t or
!! @ref fmsnetcdfunstructureddomainfile_t at the given file path name and 2D or unstructured domain.
!! @note For individual documentation on the listed routines, please see the appropriate helper module.
!! For netcdf files with a structured domain: @ref fms_netcdf_domain_io_mod.
!! For netcdf files with an unstructured domain: @ref fms_netcdf_unstructured_domain_io_mod.
!! For generic netcdf: @ref netcdf_io_mod.
!> @ingroup fms2_io_mod
interface open_file
  module procedure netcdf_file_open_wrap
  module procedure open_domain_file
  module procedure open_unstructured_domain_file
end interface open_file


!> @brief Creates a diskless netcdf or domain file. File is created in memory only via the netcdf library's
!! NC_DISKLESS creation mode option. Data will be lost upon file closing.
!!
!> @return true if successful, false otherwise
!!
!> <br>Example usage:
!!
!!              io_success = open_virtual_file(fileobj, "filename", pelist)
!!
!! Opens a virtual file through @ref fmsnetcdffile_t at an optional file path and pelist
!!
!!              io_success = open_virtual_file(fileobj, domain, "filename")
!!
!! Opens a virtual domain file through @ref fmsnetcdfdomainfile_t or
!! @ref fmsnetcdfunstructureddomainfile_t for a given 2D domain at an optional path <br>
!!
!! @note For individual documentation on the listed routines, please see the appropriate helper module: @ref blackboxio
!> @ingroup fms2_io_mod
interface open_virtual_file
  module procedure create_diskless_netcdf_file_wrap
  module procedure create_diskless_domain_file
  module procedure create_diskless_unstructured_domain_file
end interface open_virtual_file

!> @brief Close a netcdf or domain file opened with @ref open_file or
!! @ref open_virtual_file
!!
!> <br>Example usage:
!!
!!              call close_file(fileobj)
!!
!! Closes any given fileobj opened via @ref open_file or @ref open_virtual_file
!!
!! @note For individual documentation on the listed routines, please see the appropriate helper module.
!! For netcdf files with a structured domain: @ref fms_netcdf_domain_io_mod.
!! For netcdf files with an unstructured domain: @ref fms_netcdf_unstructured_domain_io_mod.
!! For generic netcdf: @ref netcdf_io_mod.
!> @ingroup fms2_io_mod
interface close_file
  module procedure netcdf_file_close_wrap
  module procedure close_domain_file
  module procedure close_unstructured_domain_file
end interface close_file

!> @brief Adds a dimension/axis to a given netcdf file object.
!!
!> <br>Example usage:
!!
!!              call register_axis(fileobj, "lon", "x")
!!
!! Adds a dimension named "lon" associated with the x axis of the 2D domain file. For unstructured
!! domains no x or y axis character is provided.
!!
!!              call register_axis(fileobj, "lon", n)
!!
!! Adds a dimension named "lon" with length n to a given netcdf file.<br>
!!
!! @note For individual documentation on the listed routines, please see the appropriate helper module.
!! For netcdf files with a structured domain: @ref fms_netcdf_domain_io_mod.
!! For netcdf files with an unstructured domain: @ref fms_netcdf_unstructured_domain_io_mod.
!! For generic netcdf: @ref netcdf_io_mod.
!> @ingroup fms2_io_mod
interface register_axis
  module procedure netcdf_add_dimension
  module procedure register_compressed_dimension
  module procedure register_domain_decomposed_dimension
  module procedure register_unstructured_dimension
end interface register_axis

!> @brief Defines a new field/variable within the given file
!> <br>Example usage:
!!
!!              call register_field(fileobj, "lon", "double", (/"lon"/) )
!!
!! Adds a double variable named "lon" to the given file, corresponding to the
!! list of dimension names (which must be previously defined in the fileobj).
!! The size of dimension name list provided is the amount of ranks for the created
!! field, scalar if list not provided.
!!
!! @note For individual documentation on the listed routines, please see the appropriate helper module.
!! For netcdf files with a structured domain: @ref fms_netcdf_domain_io_mod.
!! For netcdf files with an unstructured domain: @ref fms_netcdf_unstructured_domain_io_mod.
!! For generic netcdf: @ref netcdf_io_mod.
!> @ingroup fms2_io_mod
interface register_field
  module procedure netcdf_add_variable_wrap
  module procedure register_domain_variable
  module procedure register_unstructured_domain_variable
end interface register_field

!> @brief Similar to @ref register_field, but occupies the field with data for restarts
!> <br>Example usage:
!!
!!              call register_restart_field(fileobj, "temperature", data_ptr, (/"lon", "time"/) )
!!
!! Creates a restart variable and sets it to the values from data_ptr, corresponding to
!! the list of dimension names. Rank of data_ptr must equal the amount of corresponding dimensions.
!!
!! @note For individual documentation on the listed routines, please see the appropriate helper module.
!! For netcdf files with a structured domain: @ref fms_netcdf_domain_io_mod.
!! For netcdf files with an unstructured domain: @ref fms_netcdf_unstructured_domain_io_mod.
!! For generic netcdf: @ref netcdf_io_mod.
!> @ingroup fms2_io_mod
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
  module procedure register_restart_region_2d
  module procedure register_restart_region_3d
end interface register_restart_field

!> @brief Write data to a defined field within a file
!> <br>Example usage:
!!
!!              call write_data(fileobj, "lon", data)
!!
!! Write the value(s) in data to the field named "lon"
!!
!> @ingroup fms2_io_mod
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

!> @brief Read data from a defined field in a file
!!
!> <br>Example usage:
!!
!!              call read_data(fileobj, "lat", data)
!!
!! Read the values for the field "lat" from the file and write them onto data <br>
!!
!> @ingroup fms2_io_mod
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

!> @brief Writes all restart fields registered within a given restart file
!> <br>Example usage:
!!
!!              call write_restart(fileobj)
!!
!! Writes previously registered restart fields to the given restart file
!!
!! @note For individual documentation on the listed routines, please see the appropriate helper module.
!! For netcdf files with a structured domain: @ref fms_netcdf_domain_io_mod.
!! For netcdf files with an unstructured domain: @ref fms_netcdf_unstructured_domain_io_mod.
!! For generic netcdf: @ref netcdf_io_mod.
!> @ingroup fms2_io_mod
interface write_restart
  module procedure netcdf_save_restart_wrap
  module procedure save_domain_restart
  module procedure unstructured_write_restart
end interface write_restart

!> @brief Writes all restart fields in a given restart file to a new restart file
!> <br>Example usage:
!!
!!              call write_new_restart(fileobj, timestamp="tstring", filename="new_restartfilename")
!!
!! Creates a new restart file, with the provided timestamp and filename, out of the registered
!! restart fields in the given restart file.
!!
!! @note For individual documentation on the listed routines, please see the appropriate helper module: @ref blackboxio
!> @ingroup fms2_io_mod
interface write_new_restart
  module procedure netcdf_save_restart_wrap2
  module procedure save_domain_restart_wrap
  module procedure unstructured_write_restart_wrap
end interface write_new_restart

!> @brief Reads in restart variables from a given file
!> <br>Example usage:
!!              call read_restart(fileobj)
!! Reads registered restart variables from fileobj
!!
!! @note For individual documentation on the listed routines, please see the appropriate helper module.
!! For netcdf files with a structured domain: @ref fms_netcdf_domain_io_mod.
!! For generic netcdf: @ref netcdf_io_mod.
!> @ingroup fms2_io_mod
interface read_restart
  module procedure netcdf_restore_state
  module procedure restore_domain_state
end interface read_restart

!> @brief Read registered restarts from a new file
!! Optionally takes directory to write to, model time and filename
!> <br>Example usage:
!!              call read_new_restart(fileobj, unlimted_dimension_level)
!!
!!              call read_new_restart(fileobj, unlimited_dimension_level, directory, timestamp, filename)
!! @note For individual documentation on the listed routines, please see the appropriate helper module: @ref blackboxio
!> @ingroup fms2_io_mod
interface read_new_restart
  module procedure netcdf_restore_state_wrap
  module procedure restore_domain_state_wrap
end interface read_new_restart

!> @addtogroup fms2_io_mod
!> @{

logical, private :: fms2_io_is_initialized = .false. !< True after fms2_io_init is run
! Namelist variables
integer :: ncchksz = 64*1024  !< User defined chunksize (in bytes) argument in netcdf file
                              !! creation calls. Replaces setting the NC_CHKSZ environment variable.
character (len = 10) :: netcdf_default_format = "64bit" !< User defined netcdf file format, acceptable values
                              !! are: "64bit", "classic", "netcdf4". This can be overwritten if you specify
                              !! "nc_format" in the open_file call
integer :: header_buffer_val = 16384 !< Use defined netCDF header buffer size(in bytes) used in
                                     !! NF__ENDDEF
integer :: deflate_level = default_deflate_level !< Netcdf deflate level to use in nf90_def_var
                                                 !! (integer between 1 to 9)
logical :: shuffle = .false. !< Flag indicating whether to use the netcdf shuffle filter
namelist / fms2_io_nml / &
                      ncchksz, netcdf_default_format, header_buffer_val, deflate_level, shuffle

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
  if (header_buffer_val .le. 0) then
    call mpp_error(FATAL, "header_buffer_val in fms2_io_nml must be a positive number.")
  endif
  if (deflate_level .lt. 0 .or. deflate_level .gt. 9) then
    call mpp_error(FATAL, &
      "deflate_level in fms2_io_nml must be a positive number between 1 and 9 as it is required by NetCDF")
  endif
  call netcdf_io_init (ncchksz,header_buffer_val,netcdf_default_format, deflate_level, shuffle)
  call blackboxio_init (ncchksz)
!> Mark the fms2_io as initialized
  fms2_io_is_initialized = .true.
end subroutine fms2_io_init

end module fms2_io_mod
!> @}
! close documentation grouping
