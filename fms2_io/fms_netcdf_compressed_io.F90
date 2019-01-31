!> @file
module fms_netcdf_compressed_io_mod
use,intrinsic :: iso_fortran_env
use netcdf
use fms_io_utils_mod
use mpp_mod
use netcdf_io_mod
implicit none
private


!Module constants.
integer,parameter :: dimension_not_found = 0
integer,parameter,public :: max_num_compressed_dims = 10


!> @brief Compressed dimension.
type CompressedDimension_t
    private
    character(len=256) :: dimname !< Dimension name.
    integer,dimension(:),allocatable :: npes_corner !< Array of starting
                                                    !! indices for each rank.
    integer,dimension(:),allocatable :: npes_nelems !< Number of elements
                                                    !! associated with each
                                                    !! rank.
    integer :: nelems !< Total size of the dimension.
endtype CompressedDimension_t


!> @brief netcdf compressed file type.
type,extends(NetcdfFile_t),public :: FmsNetcdfCompressedFile_t
    type(CompressedDimension_t),dimension(:),allocatable :: compressed_dims !< "Compressed" dimension.
    integer :: n !< Number of compressed dimensions.
    type(char_linked_list),pointer :: compressed_vars => null() !< Variables that depend
                                                                !! on a compressed
                                                                !! dimension.
endtype FmsNetcdfCompressedFile_t


public :: open_compressed_file
public :: close_compressed_file
public :: register_compressed_dimension
public :: register_compressed_variable
public :: register_compressed_restart_variable_0d
public :: register_compressed_restart_variable_1d
public :: register_compressed_restart_variable_2d
public :: register_compressed_restart_variable_3d
public :: register_compressed_restart_variable_4d
public :: register_compressed_restart_variable_5d
public :: compressed_read_0d
public :: compressed_read_1d
public :: compressed_read_2d
public :: compressed_read_3d
public :: compressed_read_4d
public :: compressed_read_5d
public :: compressed_write_0d
public :: compressed_write_1d
public :: compressed_write_2d
public :: compressed_write_3d
public :: compressed_write_4d
public :: compressed_write_5d
public :: save_compressed_restart
public :: compressed_start_and_count


contains


!> @brief Get the index of a compressed dimension in a file object.
!! @return Index of the compressed dimension.
function get_compressed_dimension_index(fileobj, &
                                        dim_name) &
    result(dindex)

    !Inputs/outputs.
    class(FmsNetcdfCompressedFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: dim_name !< Dimension name.

    !Local variables.
    integer :: dindex
    integer :: i

    dindex = dimension_not_found
    do i = 1,fileobj%n
        if (string_compare(fileobj%compressed_dims(i)%dimname,dim_name)) then
            dindex = i
            return
        endif
    enddo
end function get_compressed_dimension_index


!> @brief Add a compressed dimension.
subroutine append_compressed_dimension(fileobj, &
                                       dim_name, &
                                       npes_corner, &
                                       npes_nelems)

    !Inputs/outputs.
    class(FmsNetcdfCompressedFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: dim_name !< Dimension name.
    integer,dimension(:),intent(in) :: npes_corner !< Array of starting
                                                   !! indices for each rank.
    integer,dimension(:),intent(in) :: npes_nelems !< Number of elements
                                                   !! associated with each
                                                   !! rank.

    if (get_compressed_dimension_index(fileobj,dim_name) .ne. &
        dimension_not_found) then
        call error("dimension "//trim(dim_name)//" already registered" &
                   //" to file "//trim(fileobj%path)//".")
    endif
    fileobj%n = fileobj%n + 1
    if (fileobj%n .gt. max_num_compressed_dims) then
        call error("number of compressed dimensions exceeds limit.")
    endif
    call string_copy(fileobj%compressed_dims(fileobj%n)%dimname, &
                     dim_name)
    if (size(npes_corner) .ne. size(fileobj%pelist) .or. &
        size(npes_nelems) .ne. size(fileobj%pelist)) then
        call error("incorrect size for input npes_corner or npes_nelems" &
                   //" arrays.")
    endif
    allocate(fileobj%compressed_dims(fileobj%n)%npes_corner(size(fileobj%pelist)))
    fileobj%compressed_dims(fileobj%n)%npes_corner(:) = npes_corner(:)
    allocate(fileobj%compressed_dims(fileobj%n)%npes_nelems(size(fileobj%pelist)))
    fileobj%compressed_dims(fileobj%n)%npes_nelems(:) = npes_nelems(:)
    fileobj%compressed_dims(fileobj%n)%nelems = &
        sum(fileobj%compressed_dims(fileobj%n)%npes_nelems)
end subroutine append_compressed_dimension


!> @brief Given a compressed variable, get the index of the compressed
!!        dimension.
!! @return Index of the compressed dimension.
function get_variable_compressed_dimension_index(fileobj, &
                                                 variable_name, &
                                                 broadcast) &
    result(compressed_dimension_index)

    !Inputs/outputs.
    class(FmsNetcdfCompressedFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    integer :: compressed_dimension_index

    !Local variables.
    character(len=nf90_max_name),dimension(:),allocatable :: dim_names
    integer :: i

    compressed_dimension_index = dimension_not_found
    if (fileobj%is_root) then
        call get_variable_dimension_names(fileobj, &
                                          variable_name, &
                                          dim_names, &
                                          broadcast=.false.)
        do i = 1,size(dim_names)
            if (get_compressed_dimension_index(fileobj,dim_names(i)) .ne. &
                dimension_not_found) then
                compressed_dimension_index = i
                exit
            endif
        enddo
        if (compressed_dimension_index .eq. dimension_not_found) then
            call error("no compresssed dimension found.")
        endif
        deallocate(dim_names)
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    call mpp_broadcast(compressed_dimension_index, &
                       fileobj%io_root, &
                       pelist=fileobj%pelist)
end function get_variable_compressed_dimension_index


!> @brief Add a compressed variable.
subroutine append_compressed_variable(fileobj, &
                                      variable_name, &
                                      dimensions)

    !Inputs/outputs.
    class(FmsNetcdfCompressedFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.

    !Local variables.
    logical :: found
    integer :: i

    if (present(dimensions)) then
        found = .false.
        do i = 1,size(dimensions)
            if (get_compressed_dimension_index(fileobj,dimensions(i)) .ne. &
                dimension_not_found) then
                if (found) then
                    call error("multiple compressed dimensions found.")
                else
                    found = .true.
                endif
            endif
        enddo
        if (found) then
            call append_to_list(fileobj%compressed_vars, &
                                variable_name)
        endif
    endif
end subroutine append_compressed_variable


!> @brief Open a netcdf file that will contain compressed dimensions.
!! @return .true. if open succeeds, or else .false.
function open_compressed_file(fileobj, &
                              path, &
                              mode, &
                              pelist, &
                              nc_format, &
                              is_restart) &
    result(success)

    !Inputs/outputs.
    type(FmsNetcdfCompressedFile_t),intent(inout) :: fileobj
    character(len=*),intent(in) :: path !< File path.
    character(len=*),intent(in) :: mode !< File mode.  Allowed values
                                        !! are "read", "append", "write", or
                                        !! "overwrite".
    integer,dimension(:),intent(in) :: pelist !< A list of ranks associated
                                              !! with this file.
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

    success = netcdf_file_open(fileobj, &
                               path, &
                               mode, &
                               nc_format=nc_format, &
                               pelist=pelist, &
                               is_restart=is_restart)
    if (.not. success) then
        return
    endif
    allocate(fileobj%compressed_dims(max_num_compressed_dims))
    fileobj%n = 0
end function open_compressed_file


!> @brief Close a netcdf file that contains compressed dimensions.
subroutine close_compressed_file(fileobj)

    !Inputs/outputs.
    class(FmsNetcdfCompressedFile_t),intent(inout) :: fileobj !< File object.

    !Local variables.
    integer :: i

    call netcdf_file_close(fileobj)
    do i = 1,fileobj%n
        if (allocated(fileobj%compressed_dims(i)%npes_corner)) then
            deallocate(fileobj%compressed_dims(i)%npes_corner)
        endif
        if (allocated(fileobj%compressed_dims(i)%npes_nelems)) then
            deallocate(fileobj%compressed_dims(i)%npes_nelems)
        endif
    enddo
    deallocate(fileobj%compressed_dims)
    call destroy_list(fileobj%compressed_vars)
end subroutine close_compressed_file


!> @brief Add a compressed dimension.
subroutine register_compressed_dimension(fileobj, &
                                         dim_name, &
                                         dim_size, &
                                         npes_corner, &
                                         npes_nelems)

    !Inputs/outputs.
    class(FmsNetcdfCompressedFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: dim_name !< Dimension name.
    integer,intent(in),optional :: dim_size !< Dimension size.
    integer,dimension(:),intent(in),optional :: npes_corner !< Array of starting
                                                            !! indices for each rank.
    integer,dimension(:),intent(in),optional :: npes_nelems !< Number of elements
                                                            !! associated with each
                                                            !! rank.

    !Local variables.
    integer :: dsize
    integer :: fdim_size

    if (present(npes_corner) .or. present(npes_nelems)) then
        if (.not. present(npes_corner) .or. .not. present(npes_nelems)) then
            call error("both npes_corner and npes_nelems are required in" &
                       //" order to register a compressed dimension.")
        endif
        call append_compressed_dimension(fileobj, &
                                         dim_name, &
                                         npes_corner, &
                                         npes_nelems)
        dsize = sum(npes_nelems)
        if (fileobj%is_readonly) then
            call get_dimension_size(fileobj, &
                                    dim_name, &
                                    fdim_size, &
                                    broadcast=.true.)
            if (fdim_size .ne. dsize) then
                call error("dimension "//trim(dim_name)//" does not match" &
                           //" the size of the associated compressed axis.")
            endif
        endif
    elseif (.not. fileobj%is_readonly) then
        if (.not. present(dim_size)) then
            call error("dim_size is required for non-compressed dimensions.")
        endif
        dsize = dim_size
    endif
    if (.not. fileobj%is_readonly) then
        call netcdf_add_dimension(fileobj, &
                                  dim_name, &
                                  dsize)
    endif
end subroutine register_compressed_dimension


!> @brief Add a compressed variable.
subroutine register_compressed_variable(fileobj, &
                                        variable_name, &
                                        variable_type, &
                                        dimensions)

    !Inputs/outputs.
    class(FmsNetcdfCompressedFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    character(len=*),intent(in) :: variable_type !< Variable type.  Allowed
                                                 !! values are: "int", "int64",
                                                 !! "float", or "double".
    character(len=*),dimension(:),intent(in),optional :: dimensions !< Dimension names.

    call append_compressed_variable(fileobj, &
                                    variable_name, &
                                    dimensions)
    call netcdf_add_variable(fileobj, &
                             variable_name, &
                             variable_type, &
                             dimensions)
end subroutine register_compressed_variable


!> @brief Loop through registered restart variables and write them to
!!        a netcdf file.
subroutine save_compressed_restart(fileobj, &
                                   unlim_dim_level)

    !Inputs/outputs.
    class(FmsNetcdfCompressedFile_t),intent(in) :: fileobj !< File object.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.

    !Local variables.
    integer :: i

    if (.not. fileobj%is_restart) then
        call error("file "//trim(fileobj%path)//" is not a restart file.")
    endif
    do i = 1,fileobj%num_restart_vars
        if (associated(fileobj%restart_vars(i)%data0d)) then
            call compressed_write_0d(fileobj, &
                                     fileobj%restart_vars(i)%varname, &
                                     fileobj%restart_vars(i)%data0d, &
                                     unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data1d)) then
            call compressed_write_1d(fileobj, &
                                     fileobj%restart_vars(i)%varname, &
                                     fileobj%restart_vars(i)%data1d, &
                                     unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data2d)) then
            call compressed_write_2d(fileobj, &
                                     fileobj%restart_vars(i)%varname, &
                                     fileobj%restart_vars(i)%data2d, &
                                     unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data3d)) then
            call compressed_write_3d(fileobj, &
                                     fileobj%restart_vars(i)%varname, &
                                     fileobj%restart_vars(i)%data3d, &
                                     unlim_dim_level=unlim_dim_level)
        elseif (associated(fileobj%restart_vars(i)%data4d)) then
            call compressed_write_4d(fileobj, &
                                     fileobj%restart_vars(i)%varname, &
                                     fileobj%restart_vars(i)%data4d, &
                                     unlim_dim_level=unlim_dim_level)
        else
            call error("this branch should not be reached.")
        endif
    enddo
end subroutine save_compressed_restart


!> @brief Gathers a compressed arrays size and offset for each pe.
subroutine compressed_start_and_count(fileobj, nelems, npes_start, npes_count)

  class(FmsNetcdfCompressedFile_t), intent(in) :: fileobj !< File object.
  integer, intent(in) :: nelems !< Number of elements on the current pe.
  integer, dimension(:), allocatable, intent(out) :: npes_start !< Offset for each pe.
  integer, dimension(:), allocatable, intent(out) :: npes_count !< Number of elements for
                                                                !! each pe.

  integer :: i

  allocate(npes_start(size(fileobj%pelist)))
  allocate(npes_count(size(fileobj%pelist)))
  do i = 1, size(fileobj%pelist)
    if (fileobj%pelist(i) .eq. mpp_pe()) then
      npes_count(i) = nelems
    else
      call mpp_recv(npes_count(i), fileobj%pelist(i), block=.false.)
      call mpp_send(nelems, fileobj%pelist(i))
    endif
  enddo
  call mpp_sync_self(check=EVENT_RECV)
  call mpp_sync_self(check=EVENT_SEND)
  npes_start(1) = 1
  do i = 1, size(fileobj%pelist)-1
    npes_start(i+1) = npes_start(i) + npes_count(i)
  enddo
end subroutine compressed_start_and_count


include "register_compressed_restart_variable.inc"
include "compressed_read.inc"
include "compressed_write.inc"


end module fms_netcdf_compressed_io_mod
