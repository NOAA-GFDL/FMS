!> @file

!> @brief Create an "abstract" netcdf type, which can be extended to meet
!!        our various I/O needs.
module netcdf_io_mod
use,intrinsic :: iso_fortran_env
use netcdf
use mpp_mod
use fms_io_utils_mod
implicit none
private


!Module constants.
integer,parameter :: variable_missing = -1
integer,parameter :: dimension_missing = -1
integer,parameter,public :: no_unlimited_dimension = -1
character(len=1),parameter :: missing_path = ""
integer,parameter :: missing_ncid = -1
integer,parameter :: missing_rank = -1
integer,parameter :: define_mode = 0
integer,parameter :: data_mode = 1
integer,parameter :: max_num_restart_vars = 200
integer,parameter,public :: unlimited = nf90_unlimited


!> @brief Restart variable.
type :: RestartVariable_t
    character(len=256) :: varname
    class(*),pointer :: data0d => null()
    class(*),dimension(:),pointer :: data1d => null()
    class(*),dimension(:,:),pointer :: data2d => null()
    class(*),dimension(:,:,:),pointer :: data3d => null()
    class(*),dimension(:,:,:,:),pointer :: data4d => null()
    class(*),dimension(:,:,:,:,:),pointer :: data5d => null()
endtype RestartVariable_t


!> @brief netcdf file type.
type,public :: NetcdfFile_t
    character(len=256) :: path !< File path.
    logical :: is_readonly !< Flag telling if the file is readonly.
    integer :: ncid !< Netcdf file id.
    integer,dimension(:),allocatable :: pelist !< List of ranks who will
                                               !! communicate.
    integer :: io_root !< I/O root rank of the pelist.
    logical :: is_root !< Flag telling if the current rank is the
                       !! I/O root.
    logical :: is_restart !< Flag telling if the this file is a restart
                          !! file (that has internal pointers to data).
    type(RestartVariable_t),dimension(:),allocatable :: restart_vars !< Array of registered
                                                                     !! restart variables.
    integer :: num_restart_vars !< Number of registered restart variables.
endtype NetcdfFile_t


public :: netcdf_file_open
public :: netcdf_file_close
public :: netcdf_add_dimension
public :: netcdf_add_variable
public :: netcdf_add_restart_variable
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
public :: get_variable_unlimited_dimension_index
public :: netcdf_read_data
public :: netcdf_write_data
public :: netcdf_save_restart
public :: netcdf_restore_state


interface netcdf_add_restart_variable
    module procedure netcdf_add_restart_variable_0d
    module procedure netcdf_add_restart_variable_1d
    module procedure netcdf_add_restart_variable_2d
    module procedure netcdf_add_restart_variable_3d
    module procedure netcdf_add_restart_variable_4d
    module procedure netcdf_add_restart_variable_5d
end interface netcdf_add_restart_variable


interface netcdf_read_data
    module procedure netcdf_read_data_0d
    module procedure netcdf_read_data_1d
    module procedure netcdf_read_data_2d
    module procedure netcdf_read_data_3d
    module procedure netcdf_read_data_4d
    module procedure netcdf_read_data_5d
end interface netcdf_read_data


interface netcdf_write_data
    module procedure netcdf_write_data_0d
    module procedure netcdf_write_data_1d
    module procedure netcdf_write_data_2d
    module procedure netcdf_write_data_3d
    module procedure netcdf_write_data_4d
    module procedure netcdf_write_data_5d
end interface netcdf_write_data


interface register_global_attribute
    module procedure register_global_attribute_0d
    module procedure register_global_attribute_1d
end interface register_global_attribute


interface register_variable_attribute
    module procedure register_variable_attribute_0d
    module procedure register_variable_attribute_1d
end interface register_variable_attribute


interface get_global_attribute
    module procedure get_global_attribute_0d
    module procedure get_global_attribute_1d
end interface get_global_attribute


interface get_variable_attribute
    module procedure get_variable_attribute_0d
    module procedure get_variable_attribute_1d
end interface get_variable_attribute


contains


!> @brief Check for errors returned by netcdf.
!! @internal
subroutine check_netcdf_code(err)
    integer,intent(in) :: err !< Code returned by netcdf.
    character(len=80) :: buf
    if (err .ne. nf90_noerr) then
        buf = nf90_strerror(err)
        call error(trim(buf))
    endif
end subroutine check_netcdf_code


!> @brief Switch to the correct netcdf mode.
subroutine set_netcdf_mode(ncid, &
                           mode)
    integer,intent(in) :: ncid !< Netcdf file id.
    integer,intent(in) :: mode !< Netcdf file mode.
    integer :: err
    if (mode .eq. define_mode) then
        err = nf90_redef(ncid)
        if (err .eq. nf90_eindefine .or. err .eq. nf90_eperm) then
            return
        endif
    elseif (mode .eq. data_mode) then
        err = nf90_enddef(ncid)
        if (err .eq. nf90_enotindefine .or. err .eq. nf90_eperm) then
            return
        endif
    else
        call error("mode must be either define_mode or data_mode.")
    endif
    call check_netcdf_code(err)
end subroutine set_netcdf_mode


!> @brief Get the id of a dimension from its name.
!! @return Dimension id, or dimension_missing if it doesn't exist.
!! @internal
function get_dimension_id(ncid, &
                          dimension_name, &
                          allow_failure) &
    result(dimid)
    integer,intent(in) :: ncid !< Netcdf file id.
    character(len=*),intent(in) :: dimension_name !< Dimension name.
    logical,intent(in),optional :: allow_failure !< Flag that prevents
                                                 !! crash if dimension
                                                 !! does not exist.
    integer :: dimid
    integer :: err
    err = nf90_inq_dimid(ncid, &
                         trim(dimension_name), &
                         dimid)
    if (present(allow_failure)) then
        if (allow_failure .and. err .eq. nf90_ebaddim) then
            dimid = dimension_missing
            return
        endif
    endif
    call check_netcdf_code(err)
end function get_dimension_id


!> @brief Get the id of a variable from its name.
!! @return Variable id, or variable_missing if it doesn't exist.
!! @internal
function get_variable_id(ncid, &
                         variable_name, &
                         allow_failure) &
    result(varid)
    integer,intent(in) :: ncid !< Netcdf file object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    logical,intent(in),optional :: allow_failure !< Flag that prevents
                                                 !! crash if variable does
                                                 !! not exist.
    integer :: varid
    integer :: err
    err = nf90_inq_varid(ncid, &
                         trim(variable_name), &
                         varid)
    if (present(allow_failure)) then
        if (allow_failure .and. err .eq. nf90_enotvar) then
            varid = variable_missing
            return
        endif
    endif
    call check_netcdf_code(err)
end function get_variable_id


!> @brief Add a restart variable to a NetcdfFile_t type.
subroutine add_restart_var_to_array(fileobj, &
                                    variable_name)
    class(NetcdfFile_t),intent(inout) :: fileobj !< Netcdf file object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    integer :: i
    if (.not. fileobj%is_restart) then
        call error("file "//trim(fileobj%path)//" is not a restart file.")
    endif
    do i = 1,fileobj%num_restart_vars
        if (string_compare(fileobj%restart_vars(i)%varname,variable_name,.true.)) then
            call error("variable "//trim(variable_name)//" has already" &
                       //" been added to restart file " &
                       //trim(fileobj%path)//".")
        endif
    enddo
    fileobj%num_restart_vars = fileobj%num_restart_vars + 1
    if (fileobj%num_restart_vars .gt. max_num_restart_vars) then
        call error("Number of restart variables exceeds limit.")
    endif
    call string_copy(fileobj%restart_vars(fileobj%num_restart_vars)%varname, &
                     variable_name)
end subroutine add_restart_var_to_array


!> @brief Open a netcdf file.
!! @return .true. if open succeeds, or else .false.
function netcdf_file_open(fileobj, &
                          path, &
                          mode, &
                          nc_format, &
                          pelist, &
                          is_restart) &
    result(success)

    !Inputs/outputs.
    class(NetcdfFile_t),intent(inout) :: fileobj !< File object.
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

    !Local variables.
    integer :: nc_format_param
    integer :: err
    character(len=256) :: buf
    logical :: is_res

    !Add ".res" to the file path if necessary.
    call string_copy(buf, &
                     trim(path))
    if (present(is_restart)) then
        is_res = is_restart
    else
        is_res = .false.
    endif
    if (is_res) then
        call restart_filepath_mangle(buf, &
                                     buf)
    endif

    !Check if the file exists.
    if (string_compare(mode,"read",.true.) .or. &
        string_compare(mode,"append",.true.)) then
        success = file_exists(buf)
        if (.not. success) then
            return
        endif
    else
        success = .true.
    endif

    !Store properties in the derived type.
    call string_copy(fileobj%path, &
                     trim(buf))
    if (allocated(fileobj%pelist)) then
        deallocate(fileobj%pelist)
    endif
    if (present(pelist)) then
        allocate(fileobj%pelist(size(pelist)))
        fileobj%pelist = pelist
    else
        allocate(fileobj%pelist(1))
        fileobj%pelist(1) = mpp_pe()
    endif
    fileobj%io_root = fileobj%pelist(1)
    if (mpp_pe() .eq. fileobj%io_root) then
        fileobj%is_root = .true.
    else
        fileobj%is_root = .false.
    endif
    fileobj%is_restart = is_res
    if (fileobj%is_restart) then
        allocate(fileobj%restart_vars(max_num_restart_vars))
        fileobj%num_restart_vars = 0
    endif
    if (string_compare(mode,"read",.true.)) then
        fileobj%is_readonly = .true.
    else
        fileobj%is_readonly = .false.
    endif

    !Open the file with netcdf if this rank is the I/O root.
    if (fileobj%is_root) then
        if (present(nc_format)) then
            if (string_compare(nc_format,"64bit",.true.)) then
                nc_format_param = nf90_64bit_offset
            elseif (string_compare(nc_format,"classic",.true.)) then
                nc_format_param = nf90_classic_model
            elseif (string_compare(nc_format,"netcdf4",.true.)) then
                nc_format_param = nf90_hdf5
            else
                call error("unrecognized netcdf file format " &
                           //trim(nc_format)//".")
            endif
        else
            nc_format_param = nf90_64bit_offset
        endif
        if (string_compare(mode,"read",.true.)) then
            err = nf90_open(trim(fileobj%path), &
                            nf90_nowrite, &
                            fileobj%ncid)
        elseif (string_compare(mode,"append",.true.)) then
            err = nf90_open(trim(fileobj%path), &
                            nf90_write, &
                            fileobj%ncid)
        elseif (string_compare(mode,"write",.true.)) then
            err = nf90_create(trim(fileobj%path), &
                              ior(nf90_noclobber,nc_format_param), &
                              fileobj%ncid)
        elseif (string_compare(mode,"overwrite",.true.)) then
            err = nf90_create(trim(fileobj%path), &
                              ior(nf90_clobber,nc_format_param), &
                              fileobj%ncid)
        else
            call error("unrecognized file mode "//trim(mode)//".")
        endif
        call check_netcdf_code(err)
    else
        fileobj%ncid = missing_ncid
    endif
end function netcdf_file_open


!> @brief Close a netcdf file.
subroutine netcdf_file_close(fileobj)
    class(NetcdfFile_t),intent(inout) :: fileobj !< File object.
    integer :: err
    if (fileobj%is_root) then
        err = nf90_close(fileobj%ncid)
        call check_netcdf_code(err)
    endif
    fileobj%path = missing_path
    fileobj%ncid = missing_ncid
    if (allocated(fileobj%pelist)) then
        deallocate(fileobj%pelist)
    endif
    fileobj%io_root = missing_rank
    fileobj%is_root = .false.
    if (allocated(fileobj%restart_vars)) then
        deallocate(fileobj%restart_vars)
    endif
    fileobj%is_restart = .false.
    fileobj%num_restart_vars = 0
end subroutine netcdf_file_close


!> @brief Add a dimension to a file.
subroutine netcdf_add_dimension(fileobj, &
                                dimension_name, &
                                dim_len)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: dimension_name !< Dimension name.
    integer,intent(in) :: dim_len !< Dimension length.
    integer :: err
    integer :: dimid
    if (fileobj%is_root) then
        call set_netcdf_mode(fileobj%ncid, &
                             define_mode)
        err = nf90_def_dim(fileobj%ncid, &
                           trim(dimension_name), &
                           dim_len, &
                           dimid)
        call check_netcdf_code(err)
    endif
end subroutine netcdf_add_dimension


!> @brief Add a variable to a file.
subroutine netcdf_add_variable(fileobj, &
                               variable_name, &
                               variable_type, &
                               dimensions)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    character(len=*),intent(in) :: variable_type !< Variable type.  Allowed
                                                 !! values are: "int", "int64",
                                                 !! "float", or "double".
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    integer :: err
    integer,dimension(:),allocatable :: dimids
    integer :: vtype
    integer :: varid
    integer :: i
    if (fileobj%is_root) then
        call set_netcdf_mode(fileobj%ncid, &
                             define_mode)
        if (string_compare(variable_type,"int",.true.)) then
            vtype = nf90_int
        elseif (string_compare(variable_type,"int64",.true.)) then
            vtype = nf90_int64
        elseif (string_compare(variable_type,"float",.true.)) then
            vtype = nf90_float
        elseif (string_compare(variable_type,"double",.true.)) then
            vtype = nf90_double
        else
            call error("unsupported type.")
        endif
        if (present(dimensions)) then
            allocate(dimids(size(dimensions)))
            do i = 1,size(dimids)
                dimids(i) = get_dimension_id(fileobj%ncid, &
                                             trim(dimensions(i)))
            enddo
            err = nf90_def_var(fileobj%ncid, &
                               trim(variable_name), &
                               vtype, &
                               dimids, &
                               varid)
            deallocate(dimids)
        else
            err = nf90_def_var(fileobj%ncid, &
                               trim(variable_name), &
                               vtype, &
                               varid)
        endif
        call check_netcdf_code(err)
    endif
end subroutine netcdf_add_variable


!> @brief Loop through registered restart variables and write them to
!!        a netcdf file.
subroutine netcdf_save_restart(fileobj, &
                               unlim_dim_level)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer :: i
    if (.not. fileobj%is_restart) then
        call error("file "//trim(fileobj%path)//" is not a restart file.")
    endif
    do i = 1,fileobj%num_restart_vars
        if (associated(fileobj%restart_vars(i)%data0d)) then
            call netcdf_write_data(fileobj, &
                                   fileobj%restart_vars(i)%varname, &
                                   fileobj%restart_vars(i)%data0d, &
                                   unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data1d)) then
            call netcdf_write_data(fileobj, &
                                   fileobj%restart_vars(i)%varname, &
                                   fileobj%restart_vars(i)%data1d, &
                                   unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data2d)) then
            call netcdf_write_data(fileobj, &
                                   fileobj%restart_vars(i)%varname, &
                                   fileobj%restart_vars(i)%data2d, &
                                   unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data3d)) then
            call netcdf_write_data(fileobj, &
                                   fileobj%restart_vars(i)%varname, &
                                   fileobj%restart_vars(i)%data3d, &
                                   unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data4d)) then
            call netcdf_write_data(fileobj, &
                                   fileobj%restart_vars(i)%varname, &
                                   fileobj%restart_vars(i)%data4d, &
                                   unlim_dim_level=unlim_dim_level)
        else
            call error("this branch should not be reached.")
        endif
    enddo
end subroutine netcdf_save_restart


!> @breif Loop through registered restart variables and read them from
!!        a netcdf file.
subroutine netcdf_restore_state(fileobj, &
                                unlim_dim_level)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer :: i
    if (.not. fileobj%is_restart) then
        call error("file "//trim(fileobj%path)//" is not a restart file.")
    endif
    do i = 1,fileobj%num_restart_vars
        if (associated(fileobj%restart_vars(i)%data0d)) then
            call netcdf_read_data(fileobj, &
                                  fileobj%restart_vars(i)%varname, &
                                  fileobj%restart_vars(i)%data0d, &
                                  unlim_dim_level=unlim_dim_level, &
                                  broadcast=.true.)
        elseif (associated(fileobj%restart_vars(i)%data1d)) then
            call netcdf_read_data(fileobj, &
                                  fileobj%restart_vars(i)%varname, &
                                  fileobj%restart_vars(i)%data1d, &
                                  unlim_dim_level=unlim_dim_level, &
                                  broadcast=.true.)
        elseif (associated(fileobj%restart_vars(i)%data2d)) then
            call netcdf_read_data(fileobj, &
                                  fileobj%restart_vars(i)%varname, &
                                  fileobj%restart_vars(i)%data2d, &
                                  unlim_dim_level=unlim_dim_level, &
                                  broadcast=.true.)
        elseif (associated(fileobj%restart_vars(i)%data3d)) then
            call netcdf_read_data(fileobj, &
                                  fileobj%restart_vars(i)%varname, &
                                  fileobj%restart_vars(i)%data3d, &
                                  unlim_dim_level=unlim_dim_level, &
                                  broadcast=.true.)
        elseif (associated(fileobj%restart_vars(i)%data4d)) then
            call netcdf_read_data(fileobj, &
                                  fileobj%restart_vars(i)%varname, &
                                  fileobj%restart_vars(i)%data4d, &
                                  unlim_dim_level=unlim_dim_level, &
                                  broadcast=.true.)
        else
            call error("this branch should not be reached.")
        endif
    enddo
end subroutine netcdf_restore_state


!> @brief Determine if an attribute exists.
!! @return Flag telling if the attribute exists.
!! @internal
function attribute_exists(ncid, &
                          varid, &
                          attribute_name) &
    result(att_exists)
    integer,intent(in) :: ncid !< Netcdf file id.
    integer,intent(in) :: varid !< Variable id.
    character(len=*),intent(in) :: attribute_name !< Attribute name.
    logical :: att_exists
    integer :: err
    err = nf90_inquire_attribute(ncid, &
                                 varid, &
                                 trim(attribute_name))
    if (err .eq. nf90_enotatt) then
        att_exists = .false.
    else
        call check_netcdf_code(err)
        att_exists = .true.
    endif
end function attribute_exists


!> @brief Determine if a global attribute exists.
!! @return Flag telling if a global attribute exists.
function global_att_exists(fileobj, &
                           attribute_name, &
                           broadcast) &
    result(att_exists)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: attribute_name !< Attribute name.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    logical :: att_exists
    if (fileobj%is_root) then
        att_exists = attribute_exists(fileobj%ncid, &
                                      nf90_global, &
                                      trim(attribute_name))
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(att_exists, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
end function global_att_exists


!> @brief Determine if a variable's attribute exists.
!! @return Flag telling if the variable's attribute exists.
function variable_att_exists(fileobj, &
                             variable_name, &
                             attribute_name, &
                             broadcast) &
    result(att_exists)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    character(len=*),intent(in) :: attribute_name !< Attribute name.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    logical :: att_exists
    integer :: varid
    if (fileobj%is_root) then
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name))
        att_exists = attribute_exists(fileobj%ncid, &
                                      varid, &
                                      trim(attribute_name))
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(att_exists, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
end function variable_att_exists


!> @brief Determine the number of dimensions in a file.
!! @return The number of dimensions in the file.
function get_num_dimensions(fileobj, &
                            broadcast) &
    result(ndims)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    integer :: ndims
    integer :: err
    if (fileobj%is_root) then
        err = nf90_inquire(fileobj%ncid, &
                           nDimensions=ndims)
        call check_netcdf_code(err)
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(ndims, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
end function get_num_dimensions


!> @brief Get the names of the dimensions in a file.
subroutine get_dimension_names(fileobj, &
                               names, &
                               broadcast)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),dimension(:),allocatable,intent(inout) :: names !< Names of the
                                                                     !! dimensions.
                                                                     !! This routine will
                                                                     !! reallocate this
                                                                     !! array to the
                                                                     !! correct number
                                                                     !! of elements.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    integer :: ndims
    integer :: i
    integer :: err
    if (allocated(names)) then
        deallocate(names)
    endif
    if (fileobj%is_root) then
        ndims = get_num_dimensions(fileobj, &
                                   broadcast=.false.)
        if (ndims .gt. 0) then
            allocate(names(ndims))
            names(:) = ""
            do i = 1,ndims
                err = nf90_inquire_dimension(fileobj%ncid, &
                                             i, &
                                             name=names(i))
                call check_netcdf_code(err)
            enddo
        endif
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(ndims, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
    if (ndims .gt. 0) then
        if (.not. fileobj%is_root) then
            allocate(names(ndims))
            names(:) = ""
        endif
        call mpp_broadcast(names, &
                           len(names(ndims)), &
                           fileobj%io_root, &
                           pelist=fileobj%pelist)
    endif
end subroutine get_dimension_names


!> @brief Determine if a dimension exists.
!! @return Flag telling if the dimension exists.
function dimension_exists(fileobj, &
                          dimension_name, &
                          broadcast) &
    result(dim_exists)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: dimension_name !< Dimension name.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    logical :: dim_exists
    integer :: dimid
    if (fileobj%is_root) then
        dimid = get_dimension_id(fileobj%ncid, &
                                 trim(dimension_name), &
                                 allow_failure=.true.)
        if (dimid .eq. dimension_missing) then
            dim_exists = .false.
        else
            dim_exists = .true.
        endif
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(dim_exists, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
end function dimension_exists


!> @breif Determine where or not the dimension is unlimited.
!! @return True if the dimension is unlimited, or else false.
function is_dimension_unlimited(fileobj, &
                                dimension_name, &
                                broadcast) &
    result(is_unlimited)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: dimension_name !< Dimension name.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    logical :: is_unlimited
    integer :: dimid
    integer :: err
    integer :: ulim_dimid
    if (fileobj%is_root) then
        dimid = get_dimension_id(fileobj%ncid, &
                                 trim(dimension_name))
        err = nf90_inquire(fileobj%ncid, &
                           unlimitedDimId=ulim_dimid)
        call check_netcdf_code(err)
        is_unlimited = dimid .eq. ulim_dimid
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(is_unlimited, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
end function is_dimension_unlimited


!> @brief Get the length of a dimension.
subroutine get_dimension_size(fileobj, &
                              dimension_name, &
                              dim_size, &
                              broadcast)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: dimension_name !< Dimension name.
    integer,intent(inout) :: dim_size !< Dimension size.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    integer :: dimid
    integer :: err
    if (fileobj%is_root) then
        dimid = get_dimension_id(fileobj%ncid, &
                                 trim(dimension_name))
        err = nf90_inquire_dimension(fileobj%ncid, &
                                     dimid, &
                                     len=dim_size)
        call check_netcdf_code(err)
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(dim_size, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
end subroutine get_dimension_size


!> @brief Determine the number of variables in a file.
!! @return The number of variables in the file.
function get_num_variables(fileobj, &
                           broadcast) &
    result(nvars)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    integer :: nvars
    integer :: err
    if (fileobj%is_root) then
        err = nf90_inquire(fileobj%ncid, &
                           nVariables=nvars)
        call check_netcdf_code(err)
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(nvars, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
end function get_num_variables


!> @brief Get the names of the variables in a file.
subroutine get_variable_names(fileobj, &
                              names, &
                              broadcast)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),dimension(:),allocatable,intent(inout) :: names !< Names of the
                                                                     !! variables.
                                                                     !! This routine will
                                                                     !! reallocate this
                                                                     !! array to the
                                                                     !! correct number
                                                                     !! of elements.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    integer :: nvars
    integer :: i
    integer :: err
    if (allocated(names)) then
        deallocate(names)
    endif
    if (fileobj%is_root) then
        nvars = get_num_variables(fileobj, &
                                  broadcast=.false.)
        if (nvars .gt. 0) then
            allocate(names(nvars))
            names(:) = ""
            do i = 1,nvars
                err = nf90_inquire_variable(fileobj%ncid, &
                                            i, &
                                            name=names(i))
                call check_netcdf_code(err)
            enddo
        endif
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(nvars, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
    if (nvars .gt. 0) then
        if (.not. fileobj%is_root) then
            allocate(names(nvars))
            names(:) = ""
        endif
        call mpp_broadcast(names, &
                           len(names(nvars)), &
                           fileobj%io_root, &
                           pelist=fileobj%pelist)
    endif
end subroutine get_variable_names


!> @brief Determine if a variable exists.
!! @return Flag telling if the variable exists.
function variable_exists(fileobj, &
                         variable_name, &
                         broadcast) &
    result(var_exists)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    logical :: var_exists
    integer :: varid
    if (fileobj%is_root) then
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name), &
                                allow_failure=.true.)
        if (varid .eq. variable_missing) then
            var_exists = .false.
        else
            var_exists = .true.
        endif
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(var_exists, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
end function variable_exists


!> @brief Get the number of dimensions a variable depends on.
!! @return Number of dimensions that the variable depends on.
function get_variable_num_dimensions(fileobj, &
                                     variable_name, &
                                     broadcast) &
    result(ndims)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    integer :: ndims
    integer :: varid
    integer :: err
    if (fileobj%is_root) then
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name))
        err = nf90_inquire_variable(fileobj%ncid, &
                                    varid, &
                                    ndims=ndims)
        call check_netcdf_code(err)
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(ndims, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
end function get_variable_num_dimensions


!> @brief Get the name of a variable's dimensions.
subroutine get_variable_dimension_names(fileobj, &
                                        variable_name, &
                                        dim_names, &
                                        broadcast)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    character(len=*),dimension(:),allocatable,intent(inout) :: dim_names !< Array of
                                                                         !! dimension
                                                                         !! names.
                                                                         !! This routine will
                                                                         !! reallocate this
                                                                         !! array to the
                                                                         !! correct number
                                                                         !! of elements.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    integer :: varid
    integer :: err
    integer :: ndims
    integer,dimension(nf90_max_var_dims) :: dimids
    integer :: i
    if (allocated(dim_names)) then
        deallocate(dim_names)
    endif
    if (fileobj%is_root) then
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name))
        err = nf90_inquire_variable(fileobj%ncid, &
                                    varid, &
                                    ndims=ndims, &
                                    dimids=dimids)
        call check_netcdf_code(err)
        if (ndims .gt. 0) then
            allocate(dim_names(ndims))
            dim_names(:) = ""
            do i = 1,ndims
                err = nf90_inquire_dimension(fileobj%ncid, &
                                             dimids(i), &
                                             name=dim_names(i))
                call check_netcdf_code(err)
            enddo
        endif
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(ndims, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
    if (ndims .gt. 0) then
        if (.not. fileobj%is_root) then
            allocate(dim_names(ndims))
            dim_names(:) = ""
        endif
        call mpp_broadcast(dim_names, &
                           len(dim_names(ndims)), &
                           fileobj%io_root, &
                           pelist=fileobj%pelist)
    endif
end subroutine get_variable_dimension_names


!> @brief Get the size of a variable's dimensions.
subroutine get_variable_size(fileobj, &
                             variable_name, &
                             dim_sizes, &
                             broadcast)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    integer,dimension(:),allocatable,intent(inout) :: dim_sizes !< Array of dimension
                                                                !! sizes.
                                                                !! This routine will
                                                                !! reallocate this
                                                                !! array to the
                                                                !! correct number
                                                                !! of elements.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    integer :: varid
    integer :: err
    integer :: ndims
    integer,dimension(nf90_max_var_dims) :: dimids
    integer :: i
    if (allocated(dim_sizes)) then
        deallocate(dim_sizes)
    endif
    if (fileobj%is_root) then
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name))
        err = nf90_inquire_variable(fileobj%ncid, &
                                    varid, &
                                    ndims=ndims, &
                                    dimids=dimids)
        call check_netcdf_code(err)
        if (ndims .gt. 0) then
            allocate(dim_sizes(ndims))
            do i = 1,ndims
                err = nf90_inquire_dimension(fileobj%ncid, &
                                             dimids(i), &
                                             len=dim_sizes(i))
                call check_netcdf_code(err)
            enddo
        endif
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(ndims, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
    if (ndims .gt. 0) then
        if (.not. fileobj%is_root) then
            allocate(dim_sizes(ndims))
        endif
        call mpp_broadcast(dim_sizes, &
                           ndims, &
                           fileobj%io_root, &
                           pelist=fileobj%pelist)
    endif
end subroutine get_variable_size


!> @brief Get the index of a variable's unlimited dimensions.
!! @return The index of the unlimited dimension, or else
!!         no_unlimited_dimension.
function get_variable_unlimited_dimension_index(fileobj, &
                                                variable_name, &
                                                broadcast) &
    result(unlim_dim_index)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    integer :: unlim_dim_index
    character(len=nf90_max_name),dimension(:),allocatable :: dim_names
    integer :: i
    unlim_dim_index = no_unlimited_dimension
    if (fileobj%is_root) then
        call get_variable_dimension_names(fileobj, &
                                          variable_name, &
                                          dim_names, &
                                          broadcast=.false.)
        do i = 1,size(dim_names)
            if (is_dimension_unlimited(fileobj,dim_names(i),.false.)) then
                unlim_dim_index = i
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
    call mpp_broadcast(unlim_dim_index, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
end function get_variable_unlimited_dimension_index


include "netcdf_add_restart_variable.inc"
include "netcdf_read_data.inc"
include "netcdf_write_data.inc"
include "register_global_attribute.inc"
include "register_variable_attribute.inc"
include "get_global_attribute.inc"
include "get_variable_attribute.inc"


end module netcdf_io_mod
