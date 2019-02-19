!> @file
module fms_netcdf_io_mod
use,intrinsic :: iso_fortran_env
use fms_io_utils_mod
use netcdf_io_mod
implicit none
private


!> @brief netcdf file type.
type,extends(NetcdfFile_t),public :: FmsNetcdfFile_t
endtype FmsNetcdfFile_t


public :: open_netcdf_file
public :: close_netcdf_file
public :: register_dimension
public :: register_variable
public :: register_restart_variable_0d
public :: register_restart_variable_1d
public :: register_restart_variable_2d
public :: register_restart_variable_3d
public :: register_restart_variable_4d
public :: register_restart_variable_5d
public :: write_data_0d
public :: write_data_1d
public :: write_data_2d
public :: write_data_3d
public :: write_data_4d
public :: write_data_5d
public :: read_data_0d
public :: read_data_1d
public :: read_data_2d
public :: read_data_3d
public :: read_data_4d
public :: read_data_5d
public :: save_restart
public :: restore_state


contains


!> @brief Wrapper for netcdf file open.
!! @return .true. if open succeeds, or else .false.
function open_netcdf_file(fileobj, &
                          path, &
                          mode, &
                          nc_format, &
                          pelist, &
                          is_restart) &
    result(success)
    type(FmsNetcdfFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: path !< File path.
    character(len=*),intent(in) :: mode !< File mode.  Allowed values are:
                                        !! "read", "append", "write", or
                                        !! "overwrite".
    character(len=*),intent(in),optional :: nc_format !< Netcdf format that
                                                      !! new files are written
                                                      !! as.  Allowed values
                                                      !! are: "64bit", "classic",
                                                      !! or "netcdf4". Defaults to
                                                      !! "64bit".
    integer,dimension(:),intent(in),optional :: pelist !< List of ranks associated
                                                       !! with this file.  If not
                                                       !! provided, only the current
                                                       !! rank will be able to
                                                       !! act on the file.
    logical,intent(in),optional :: is_restart !< Flag telling if this file
                                              !! is a restart file.  Defaults
                                              !! to false.
    logical :: success
    success = netcdf_file_open(fileobj, &
                               path, &
                               mode, &
                               nc_format, &
                               pelist, &
                               is_restart)
end function open_netcdf_file


!> @brief Wrapper for netcdf file close.
subroutine close_netcdf_file(fileobj)
    type(FmsNetcdfFile_t),intent(inout) :: fileobj !< File object.
    call netcdf_file_close(fileobj)
end subroutine close_netcdf_file


!> @brief Wrappr for adding a dimension to a file.
subroutine register_dimension(fileobj, &
                              dimension_name, &
                              dim_len)
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: dimension_name !< Dimension name.
    integer,intent(in) :: dim_len !< Dimension length.
    call netcdf_add_dimension(fileobj, &
                              dimension_name, &
                              dim_len)
end subroutine register_dimension


!> @brief Wrapper for adding a variable to a file.
subroutine register_variable(fileobj, &
                             variable_name, &
                             variable_type, &
                             dimensions)
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    character(len=*),intent(in) :: variable_type !< Variable type.  Allowed
                                                 !! values are: "int", "int64",
                                                 !! "float", or "double".
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    call netcdf_add_variable(fileobj, &
                             variable_name, &
                             variable_type, &
                             dimensions)
end subroutine register_variable


!> @breif Loop through registered restart variables and write them to
!!        a netcdf file.
subroutine save_restart(fileobj, &
                        unlim_dim_level)
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    call netcdf_save_restart(fileobj, &
                             unlim_dim_level)
end subroutine save_restart


!> @breif Loop through registered restart variables and read them from
!!        a netcdf file.
subroutine restore_state(fileobj, &
                         unlim_dim_level)
    type(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    call netcdf_restore_state(fileobj, &
                              unlim_dim_level)
end subroutine restore_state


include "register_restart_variable.inc"
include "write_data.inc"
include "read_data.inc"


end module fms_netcdf_io_mod
