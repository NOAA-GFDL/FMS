!> @file
module fms_netcdf_domain_io_mod
use,intrinsic :: iso_fortran_env
use netcdf
use mpp_mod
use mpp_domains_mod
use fms_io_utils_mod
use netcdf_io_mod
implicit none
private


!Module constants.
integer,parameter :: no_domain_decomposed_dimension = 0
integer,parameter :: max_num_domain_decomposed_vars = 200
integer,parameter :: variable_not_found = 0
integer,parameter :: default_domain_position = center
character(len=16),parameter :: domain_pos_att = "domain_position"
character(len=16),parameter :: domain_axis_att_name = "domain_axis"
character(len=16),parameter :: x = "x"
character(len=16),parameter :: y = "y"


!> @brief Domain variable.
type :: DomainVariable_t
    private
    character(len=256) :: varname !< Variable name.
    integer :: pos !< Domain position.
endtype DomainVariable_t


!> @brief netcdf domain file type.
type,extends(NetcdfFile_t),public :: FmsNetcdfDomainFile_t
    private
    type(domain2d) :: domain !< Two-dimensional domain.
    type(char_linked_list),pointer :: xdims => null() !< Dimensions associated
                                                      !! with the "x" axis
                                                      !! of a 2d domain.
    type(char_linked_list),pointer :: ydims => null() !< Dimensions associated
                                                      !! with the "y" axis
                                                      !! of a 2d domain.
    type(DomainVariable_t),dimension(:),allocatable :: domain_decomposed_vars !< Variables that depend
                                                                              !! on dimensions associated
                                                                              !! with a 2d domain.
    integer :: n !< Number of domain decomposed variables.
endtype FmsNetcdfDomainFile_t


public :: open_domain_file
public :: close_domain_file
public :: register_non_domain_decomposed_dimension
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


contains


!> @brief Given a variable, get the index of the "x" or "y" domain decomposed
!!        dimension.
!! @return Index of the domain decomposed dimension or else
!!         no_domain_decomposed_dimension.
function get_domain_decomposed_dimension_index(fileobj, &
                                               variable_name, &
                                               xory, &
                                               broadcast) &
    result(dim_index)

    !Inputs/outputs.
    type(FmsNetcdfDomainFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    character(len=*),intent(in) :: xory !< String telling which dimension to
                                        !! look for.  Valid values are "x"
                                        !! or "y".
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the index will be
                                             !! broadcast to non "I/O root"
                                             !! ranks.  The broadcast will
                                             !! be done by default.
    integer :: dim_index !< Index of the decomposed dimension.

    !Local variables.
    character(len=nf90_max_name),dimension(:),allocatable :: dim_names
    type(char_linked_list),pointer :: p
    integer :: i

    dim_index = no_domain_decomposed_dimension
    if (fileobj%is_root) then
        call get_variable_dimension_names(fileobj, &
                                          variable_name, &
                                          dim_names, &
                                          broadcast=.false.)
        if (string_compare(xory,x,.true.)) then
            p => fileobj%xdims
        elseif (string_compare(xory,y,.true.)) then
            p => fileobj%ydims
        else
            call error("unrecognized xory flag value.")
        endif
        do i = 1,size(dim_names)
            if (is_in_list(p,dim_names(i),.false.)) then
                dim_index = i
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
    call mpp_broadcast(dim_index, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
end function get_domain_decomposed_dimension_index


!> @brief Determine if a variable is "domain decomposed."
!! @return Flag telling if the variable is "domain decomposed."
function is_variable_domain_decomposed(fileobj, &
                                       variable_name, &
                                       broadcast, &
                                       xindex, &
                                       yindex) &
    result(is_decomposed)

    !Inputs/outputs.
    type(FmsNetcdfDomainFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the index will be
                                             !! broadcast to non "I/O root"
                                             !! ranks.  The broadcast will
                                             !! be done by default.
    integer,intent(out),optional :: xindex !< The index of the domain
                                           !! x dimension.
    integer,intent(out),optional :: yindex !< The index of the domain
                                           !! y dimension.
    logical :: is_decomposed !< Flag telling if the variable is
                             !! "domain decomposed."

    !Local variables.
    integer,dimension(2) :: indices

    indices(1) = get_domain_decomposed_dimension_index(fileobj, &
                                                       variable_name, &
                                                       x, &
                                                       broadcast)
    if (present(xindex)) then
        xindex = indices(1)
    endif
    indices(2) = get_domain_decomposed_dimension_index(fileobj, &
                                                       variable_name, &
                                                       y, &
                                                       broadcast)
    if (present(yindex)) then
        yindex = indices(2)
    endif
    is_decomposed = (indices(1) .ne. no_domain_decomposed_dimension) .and. &
                    (indices(2) .ne. no_domain_decomposed_dimension)
end function is_variable_domain_decomposed


!> @brief Get the index of a domain decomposed variable.
!! @return Index of domain decomposed variable.
function get_domain_decomposed_variable_index(fileobj, &
                                              variable_name) &
    result(vindex)

    !Inputs/outputs.
    type(FmsNetcdfDomainFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.

    !Local variables.
    integer :: vindex
    integer :: i

    vindex = variable_not_found
    do i = 1,fileobj%n
        if (string_compare(fileobj%domain_decomposed_vars(i)%varname, &
                           variable_name)) then
            vindex = i
            return
        endif
    enddo
end function get_domain_decomposed_variable_index


!> @brief Add a domain decomposed variable to the list.
subroutine append_domain_decomposed_variable(fileobj, &
                                             variable_name, &
                                             domain_position)

    !Inputs/outputs.
    type(FmsNetcdfDomainFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    integer,intent(in),optional :: domain_position

    !Local variables.
    integer :: var_index

    if (is_variable_domain_decomposed(fileobj,variable_name,.true.)) then
        var_index = get_domain_decomposed_variable_index(fileobj, &
                                                         variable_name)
        if (var_index .ne. variable_not_found) then
            call error("variable "//trim(variable_name)//" already" &
                       //" registered.")
        endif
        fileobj%n = fileobj%n + 1
        if (fileobj%n .gt. max_num_domain_decomposed_vars) then
            call error("number of domain decomposed variables exceeds" &
                       //" limit.")
        endif
        call string_copy(fileobj%domain_decomposed_vars(fileobj%n)%varname, &
                         variable_name)
        if (present(domain_position)) then
            fileobj%domain_decomposed_vars(fileobj%n)%pos = domain_position
        else
            fileobj%domain_decomposed_vars(fileobj%n)%pos = default_domain_position
        endif
    endif
end subroutine append_domain_decomposed_variable


!> @brief Given a domain decomposed variable, get its domain position.
!! @return Position of the domain decomposed variable.
function get_domain_decomposed_variable_position(fileobj, &
                                                 variable_name) &
    result(dpos)

    !Inputs/outputs.
    type(FmsNetcdfDomainFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    integer :: dpos

    dpos = get_domain_decomposed_variable_index(fileobj, &
                                                variable_name)
    if (dpos .ne. variable_not_found) then
        dpos = fileobj%domain_decomposed_vars(dpos)%pos
    else
        dpos = default_domain_position
    endif
end function get_domain_decomposed_variable_position


!> @brief Open a domain netcdf file.
!! @return Flag telling if the open completed successfully.
function open_domain_file(fileobj, &
                          path, &
                          mode, &
                          domain, &
                          nc_format, &
                          is_restart) &
    result(success)

    !Inputs/outputs.
    type(FmsNetcdfDomainFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: path !< File path.
    character(len=*),intent(in) :: mode !< File mode.  Allowed values
                                        !! are "read", "append", "write", or
                                        !! "overwrite".
    type(domain2d),intent(in) :: domain !< Two-dimensional domain.
    character(len=*),intent(in),optional :: nc_format !< Netcdf format that
                                                      !! new files are written
                                                      !! as.  Allowed values
                                                      !! are: "64bit", "classic",
                                                      !! or "netcdf4". Defaults to
                                                      !! "64bit".
    logical,intent(in),optional :: is_restart !< Flag telling if this file
                                              !! is a restart file.  Defaults
                                              !! to false.
    logical :: success

    !Local variables.
    type(domain2d),pointer :: io_domain
    integer :: pelist_size
    integer,dimension(:),allocatable :: pelist
    integer,dimension(2) :: io_layout
    character(len=256) :: buf
    integer,dimension(1) :: tile_id

    !Make sure the input domain has an I/O domain and get its pelist.
    io_domain => mpp_get_io_domain(domain)
    if (.not. associated(io_domain)) then
        call error("input domain does not have an io_domain.")
    endif
    pelist_size = mpp_get_domain_npes(io_domain)
    allocate(pelist(pelist_size))
    call mpp_get_pelist(io_domain, &
                        pelist)
    io_layout = mpp_get_io_domain_layout(domain)

    !Add the domain tile id to the file name (if necessary).
    call string_copy(buf, &
                     path)
    if (mpp_get_ntile_count(domain) .gt. 1) then
        tile_id = mpp_get_tile_id(domain)
        call domain_tile_filepath_mangle(buf, &
                                         buf, &
                                         tile_id(1))
    endif

    success = .false.
    if (string_compare(mode,"read",.true.) .or. &
        string_compare(mode,"append",.true.)) then

        !Only for reading: attempt to open non-distributed files.
        success = netcdf_file_open(fileobj, &
                                   buf, &
                                   mode, &
                                   nc_format, &
                                   pelist, &
                                   is_restart)
    endif
    if (.not. success) then

        !Add the domain tile id to the file name (if necessary).
        if (io_layout(1)*io_layout(2) .gt. 1) then
            tile_id = mpp_get_tile_id(io_domain)
            call io_domain_tile_filepath_mangle(buf, &
                                                buf, &
                                                tile_id(1))
        endif

        !Open distributed files.
        success = netcdf_file_open(fileobj, &
                                   buf, &
                                   mode, &
                                   nc_format, &
                                   pelist, &
                                   is_restart)
    endif
    deallocate(pelist)

    !This branch should only be entered if the file attempting to be read
    !does not exist.
    if (.not. success) then
         return
    endif

    !Store/initialize necessary properties.
    fileobj%domain = domain
    allocate(fileobj%domain_decomposed_vars(max_num_domain_decomposed_vars))
    fileobj%n = 0
end function open_domain_file


!> @brief Close a domain netcdf file.
subroutine close_domain_file(fileobj)

    !Inputs/outputs.
    type(FmsNetcdfDomainFile_t),intent(inout) :: fileobj !< File object.

    call netcdf_file_close(fileobj)
    call destroy_list(fileobj%xdims)
    call destroy_list(fileobj%ydims)
    deallocate(fileobj%domain_decomposed_vars)
    fileobj%n = 0
end subroutine close_domain_file


!> @brief Wrappr for adding a non-domain decomposed dimension to a file.
subroutine register_non_domain_decomposed_dimension(fileobj, &
                                                    dimension_name, &
                                                    dim_len)

    !Inputs/outputs.
    type(FmsNetcdfDomainFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: dimension_name !< Dimension name.
    integer,intent(in) :: dim_len !< Dimension length.

    call netcdf_add_dimension(fileobj, &
                              dimension_name, &
                              dim_len)
end subroutine register_non_domain_decomposed_dimension


!> @brief Add a dimension to a file associated with a two-dimensional domain.
subroutine register_domain_decomposed_dimension(fileobj, &
                                                dim_name, &
                                                xory, &
                                                domain_position)

    !Inputs/outputs
    type(FmsNetcdfDomainFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: dim_name !< Dimension name.
    character(len=*),intent(in) :: xory !< Flag telling if the dimension
                                        !! is associated with the "x" or "y"
                                        !! axis of the 2d domain.  Allowed
                                        !! values are "x" or "y".
    integer,intent(in),optional :: domain_position !< Domain position.

    !Local variables.
    integer :: dpos
    type(domain2d),pointer :: io_domain
    integer :: domain_size
    integer :: dim_size

    if (present(domain_position)) then
        dpos = domain_position
    else
        dpos = default_domain_position
    endif
    io_domain => mpp_get_io_domain(fileobj%domain)
    if (string_compare(xory,x,.true.)) then
        call append_to_list(fileobj%xdims, &
                            dim_name)
        call mpp_get_global_domain(io_domain, &
                                   xsize=domain_size, &
                                   position=dpos)
    elseif (string_compare(xory,y,.true.)) then
        call append_to_list(fileobj%ydims, &
                            dim_name)
        call mpp_get_global_domain(io_domain, &
                                   ysize=domain_size, &
                                   position=dpos)
    else
        call error("unrecognized xory flag value.")
    endif
    if (fileobj%is_readonly) then
        call get_dimension_size(fileobj, &
                                dim_name, &
                                dim_size, &
                                broadcast=.true.)
        if (dim_size .ne. domain_size) then
            call error("dimension "//trim(dim_name)//" does not match" &
                       //" the size of the associated domain "//trim(xory) &
                       //" axis.")
        endif
    else
        call netcdf_add_dimension(fileobj, &
                                  dim_name, &
                                  domain_size)
    endif
end subroutine register_domain_decomposed_dimension


!> @brief Add a domain decomposed variable.
subroutine register_domain_variable(fileobj, &
                                    variable_name, &
                                    variable_type, &
                                    dimensions, &
                                    domain_position)

    !Inputs/outputs.
    type(FmsNetcdfDomainFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    character(len=*),intent(in) :: variable_type !< Variable type.  Allowed
                                                 !! values are: "int", "int64",
                                                 !! "float", or "double".
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.
    integer,intent(in),optional :: domain_position !< Domain position.

    if (.not. fileobj%is_readonly) then
        call netcdf_add_variable(fileobj, &
                                 variable_name, &
                                 variable_type, &
                                 dimensions)
    endif
    call append_domain_decomposed_variable(fileobj, &
                                           variable_name, &
                                           domain_position)
end subroutine register_domain_variable


!> @brief Loop through registered restart variables and write them to
!!        a netcdf file.
subroutine save_domain_restart(fileobj, &
                               unlim_dim_level)

    !Inputs/outputs.
    type(FmsNetcdfDomainFile_t),intent(in) :: fileobj !< File object.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.

    !Local variables.
    integer :: i

    if (.not. fileobj%is_restart) then
        call error("file "//trim(fileobj%path)//" is not a restart file.")
    endif
    do i = 1,fileobj%num_restart_vars
        if (associated(fileobj%restart_vars(i)%data0d)) then
            call domain_write_0d(fileobj, &
                                 fileobj%restart_vars(i)%varname, &
                                 fileobj%restart_vars(i)%data0d, &
                                 unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data1d)) then
            call domain_write_1d(fileobj, &
                                 fileobj%restart_vars(i)%varname, &
                                 fileobj%restart_vars(i)%data1d, &
                                 unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data2d)) then
            call domain_write_2d(fileobj, &
                                 fileobj%restart_vars(i)%varname, &
                                 fileobj%restart_vars(i)%data2d, &
                                 unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data3d)) then
            call domain_write_3d(fileobj, &
                                 fileobj%restart_vars(i)%varname, &
                                 fileobj%restart_vars(i)%data3d, &
                                 unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data4d)) then
            call domain_write_4d(fileobj, &
                                 fileobj%restart_vars(i)%varname, &
                                 fileobj%restart_vars(i)%data4d, &
                                 unlim_dim_level=unlim_dim_level)
        else
            call error("this branch should not be reached.")
        endif
    enddo
end subroutine save_domain_restart


!> @brief Loop through registered restart variables and read them from
!!        a netcdf file.
subroutine restore_domain_state(fileobj, &
                                unlim_dim_level)

    !Inputs/outputs.
    type(FmsNetcdfDomainFile_t),intent(in) :: fileobj !< File object.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.

    !Local variables.
    integer :: i

    if (.not. fileobj%is_restart) then
        call error("file "//trim(fileobj%path)//" is not a restart file.")
    endif
    do i = 1,fileobj%num_restart_vars
        if (associated(fileobj%restart_vars(i)%data0d)) then
            call domain_read_0d(fileobj, &
                                fileobj%restart_vars(i)%varname, &
                                fileobj%restart_vars(i)%data0d, &
                                unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data1d)) then
            call domain_read_1d(fileobj, &
                                fileobj%restart_vars(i)%varname, &
                                fileobj%restart_vars(i)%data1d, &
                                unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data2d)) then
            call domain_read_2d(fileobj, &
                                fileobj%restart_vars(i)%varname, &
                                fileobj%restart_vars(i)%data2d, &
                                unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data3d)) then
            call domain_read_3d(fileobj, &
                                fileobj%restart_vars(i)%varname, &
                                fileobj%restart_vars(i)%data3d, &
                                unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data4d)) then
            call domain_read_4d(fileobj, &
                                fileobj%restart_vars(i)%varname, &
                                fileobj%restart_vars(i)%data4d, &
                                unlim_dim_level=unlim_dim_level)
        else
            call error("this branch should not be reached.")
        endif
    enddo
end subroutine restore_domain_state


include "register_domain_restart_variable.inc"
include "domain_read.inc"
include "domain_write.inc"


end module fms_netcdf_domain_io_mod
