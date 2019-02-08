module fms2_io_mod
use netcdf_io_mod, only: unlimited, &
                         NetcdfFile_t, &
                         global_att_exists, &
                         variable_att_exists, &
                         register_global_attribute, &
                         register_variable_attribute, &
                         get_global_attribute, &
                         get_variable_attribute, &
                         get_num_dimensions, &
                         get_dimension_names, &
                         dimension_exists, &
                         is_dimension_unlimited, &
                         get_dimension_size, &
                         get_num_variables, &
                         get_variable_names, &
                         variable_exists, &
                         get_variable_num_dimensions, &
                         get_variable_dimension_names, &
                         get_variable_size, &
                         get_variable_unlimited_dimension_index, &
                         Valid_t, get_valid, is_valid, &
                         get_unlimited_dimension_name
use fms_netcdf_io_mod, only: FmsNetcdfFile_t, &
                             open_netcdf_file, &
                             close_netcdf_file, &
                             register_dimension, &
                             register_variable, &
                             register_restart_variable_0d, &
                             register_restart_variable_1d, &
                             register_restart_variable_2d, &
                             register_restart_variable_3d, &
                             register_restart_variable_4d, &
                             register_restart_variable_5d, &
                             write_data_0d, &
                             write_data_1d, &
                             write_data_2d, &
                             write_data_3d, &
                             write_data_4d, &
                             write_data_5d, &
                             read_data_0d, &
                             read_data_1d, &
                             read_data_2d, &
                             read_data_3d, &
                             read_data_4d, &
                             read_data_5d, &
                             save_restart, &
                             restore_state
use fms_netcdf_domain_io_mod, only: FmsNetcdfDomainFile_t, &
                                    open_domain_file, &
                                    close_domain_file, &
                                    register_non_domain_decomposed_dimension, &
                                    register_domain_decomposed_dimension, &
                                    register_domain_variable, &
                                    register_domain_restart_variable_0d, &
                                    register_domain_restart_variable_1d, &
                                    register_domain_restart_variable_2d, &
                                    register_domain_restart_variable_3d, &
                                    register_domain_restart_variable_4d, &
                                    register_domain_restart_variable_5d, &
                                    domain_read_0d, &
                                    domain_read_1d, &
                                    domain_read_2d, &
                                    domain_read_3d, &
                                    domain_read_4d, &
                                    domain_read_5d, &
                                    domain_write_0d, &
                                    domain_write_1d, &
                                    domain_write_2d, &
                                    domain_write_3d, &
                                    domain_write_4d, &
                                    domain_write_5d, &
                                    save_domain_restart, &
                                    restore_domain_state, &
                                    add_domain_decomposition_attribute, &
                                    get_compute_domain_dimension_indices
use fms_netcdf_compressed_io_mod, only: FmsNetcdfCompressedFile_t, &
                                        open_compressed_file, &
                                        close_compressed_file, &
                                        register_compressed_dimension, &
                                        register_compressed_variable, &
                                        register_compressed_restart_variable_0d, &
                                        register_compressed_restart_variable_1d, &
                                        register_compressed_restart_variable_2d, &
                                        register_compressed_restart_variable_3d, &
                                        register_compressed_restart_variable_4d, &
                                        register_compressed_restart_variable_5d, &
                                        compressed_read_0d, &
                                        compressed_read_1d, &
                                        compressed_read_2d, &
                                        compressed_read_3d, &
                                        compressed_read_4d, &
                                        compressed_read_5d, &
                                        compressed_write_0d, &
                                        compressed_write_1d, &
                                        compressed_write_2d, &
                                        compressed_write_3d, &
                                        compressed_write_4d, &
                                        compressed_write_5d, &
                                        save_compressed_restart, &
                                        compressed_start_and_count
use fms_netcdf_unstructured_domain_io_mod, only: FmsNetcdfUnstructuredDomainFile_t, &
                                                 open_unstructured_domain_file, &
                                                 register_unstructured_dimension
implicit none
private


public :: unlimited
public :: NetcdfFile_t
public :: FmsNetcdfFile_t
public :: FmsNetcdfDomainFile_t
public :: FmsNetcdfCompressedFile_t
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
public :: compressed_start_and_count
public :: add_domain_decomposition_attribute
public :: get_compute_domain_dimension_indices
public :: Valid_t
public :: get_valid
public :: is_valid
public :: get_unlimited_dimension_name

interface open_file
    module procedure open_netcdf_file
    module procedure open_domain_file
    module procedure open_compressed_file
    module procedure open_unstructured_domain_file
end interface open_file


interface close_file
    module procedure close_netcdf_file
    module procedure close_domain_file
    module procedure close_compressed_file
end interface close_file


interface register_axis
    module procedure register_dimension
    module procedure register_non_domain_decomposed_dimension
    module procedure register_domain_decomposed_dimension
    module procedure register_compressed_dimension
    module procedure register_unstructured_dimension
end interface register_axis


interface register_field
    module procedure register_variable
    module procedure register_domain_variable
    module procedure register_compressed_variable
end interface register_field


interface register_restart_field
    module procedure register_restart_variable_0d
    module procedure register_restart_variable_1d
    module procedure register_restart_variable_2d
    module procedure register_restart_variable_3d
    module procedure register_restart_variable_4d
    module procedure register_restart_variable_5d
    module procedure register_domain_restart_variable_0d
    module procedure register_domain_restart_variable_1d
    module procedure register_domain_restart_variable_2d
    module procedure register_domain_restart_variable_3d
    module procedure register_domain_restart_variable_4d
    module procedure register_domain_restart_variable_5d
    module procedure register_compressed_restart_variable_0d
    module procedure register_compressed_restart_variable_1d
    module procedure register_compressed_restart_variable_2d
    module procedure register_compressed_restart_variable_3d
    module procedure register_compressed_restart_variable_4d
    module procedure register_compressed_restart_variable_5d
end interface register_restart_field


interface write_data
    module procedure write_data_0d
    module procedure write_data_1d
    module procedure write_data_2d
    module procedure write_data_3d
    module procedure write_data_4d
    module procedure write_data_5d
    module procedure domain_write_0d
    module procedure domain_write_1d
    module procedure domain_write_2d
    module procedure domain_write_3d
    module procedure domain_write_4d
    module procedure domain_write_5d
    module procedure compressed_write_0d
    module procedure compressed_write_1d
    module procedure compressed_write_2d
    module procedure compressed_write_3d
    module procedure compressed_write_4d
    module procedure compressed_write_5d
end interface write_data


interface read_data
    module procedure read_data_0d
    module procedure read_data_1d
    module procedure read_data_2d
    module procedure read_data_3d
    module procedure read_data_4d
    module procedure read_data_5d
    module procedure domain_read_0d
    module procedure domain_read_1d
    module procedure domain_read_2d
    module procedure domain_read_3d
    module procedure domain_read_4d
    module procedure domain_read_5d
    module procedure compressed_read_0d
    module procedure compressed_read_1d
    module procedure compressed_read_2d
    module procedure compressed_read_3d
    module procedure compressed_read_4d
    module procedure compressed_read_5d
end interface read_data


interface write_restart
    module procedure save_restart
    module procedure save_domain_restart
    module procedure save_compressed_restart
end interface write_restart


interface read_restart
    module procedure restore_state
    module procedure restore_domain_state
end interface read_restart


end module fms2_io_mod
