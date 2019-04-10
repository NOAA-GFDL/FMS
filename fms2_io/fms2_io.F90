!> @file

!> @brief Public API for fms I/O.
module fms2_io_mod
use fms_io_utils_mod
use netcdf_io_mod
use fms_netcdf_domain_io_mod
use fms_netcdf_unstructured_domain_io_mod
implicit none
private


public :: unlimited
public :: FmsNetcdfFile_t
public :: FmsNetcdfDomainFile_t
public :: FmsNetcdfUnstructuredDomainFile_t
public :: open_file
public :: close_file
public :: register_axis
public :: register_field
public :: register_restart_field
public :: write_data
public :: read_data
public :: write_restart
public :: read_restart
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
public :: Valid_t
public :: get_valid
public :: is_valid
public :: get_unlimited_dimension_name
public :: file_exists
public :: compressed_start_and_count
public :: is_registered_to_restart


interface open_file
  module procedure netcdf_file_open
  module procedure open_domain_file
  module procedure open_unstructured_domain_file
end interface open_file


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


interface read_restart
  module procedure netcdf_restore_state
  module procedure restore_domain_state
end interface read_restart


end module fms2_io_mod
