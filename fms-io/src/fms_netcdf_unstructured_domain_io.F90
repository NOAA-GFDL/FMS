!> @file
module fms_netcdf_unstructured_domain_io_mod
use,intrinsic :: iso_fortran_env
use netcdf
use mpp_domains_mod
use fms_io_utils_mod
use netcdf_io_mod
use fms_netcdf_compressed_io_mod
implicit none
private


!> @brief netcdf unstructured domain file type.
type,extends(FmsNetcdfCompressedFile_t),public :: FmsNetcdfUnstructuredDomainFile_t
    private
    type(domainug) :: domain
endtype FmsNetcdfUnstructuredDomainFile_t


public :: open_unstructured_domain_file
public :: register_unstructured_dimension

contains


!> @brief Open a netcdf file that is associated with an unstructured domain.
!! @return Flag telling if the open completed successfully.
function open_unstructured_domain_file(fileobj, &
                                       path, &
                                       mode, &
                                       domain, &
                                       nc_format, &
                                       is_restart) &
    result(success)

    !Inputs/outputs.
    type(FmsNetcdfUnstructuredDomainFile_t),intent(inout) :: fileobj
    character(len=*),intent(in) :: path !< File path.
    character(len=*),intent(in) :: mode !< File mode.  Allowed values
                                        !! are "read", "append", "write", or
                                        !! "overwrite".
    type(domainug),intent(in) :: domain !< Unstructured domain.
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
    type(domainug),pointer :: io_domain
    integer :: pelist_size
    integer,dimension(:),allocatable :: pelist
    character(len=256) :: buf
    integer :: tile_id

    !Get the input domain's I/O domain pelist.
    io_domain => mpp_get_ug_io_domain(domain)
    if (.not. associated(io_domain)) then
        call error("input domain does not have an io_domain.")
    endif
    pelist_size = mpp_get_ug_domain_npes(io_domain)
    allocate(pelist(pelist_size))
    call mpp_get_ug_domain_pelist(io_domain, &
                                  pelist)

    !Add the domain tile id to the file name (if necessary).
    call string_copy(buf, &
                     path)
    if (mpp_get_UG_domain_ntiles(domain) .gt. 1) then
        tile_id = mpp_get_ug_domain_tile_id(domain)
        call domain_tile_filepath_mangle(buf, &
                                         path, &
                                         tile_id)
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
        if (mpp_get_io_domain_ug_layout(domain) .gt. 1) then
            tile_id = mpp_get_ug_domain_tile_id(io_domain)
            call io_domain_tile_filepath_mangle(buf, &
                                                buf, &
                                                tile_id)
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
    allocate(fileobj%compressed_dims(max_num_compressed_dims))
    fileobj%n = 0
end function open_unstructured_domain_file


!> @brief Add an unstructured dimension.
subroutine register_unstructured_dimension(fileobj, &
                                           dim_name, &
                                           is_unstructured, &
                                           dim_size)
    type(FmsNetcdfUnstructuredDomainFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: dim_name !< Dimension name.
    logical,intent(in) :: is_unstructured !< Flag telling if the
                                          !! dimension is unstructured.
    integer,intent(in),optional :: dim_size !< Dimension size.
    type(domainug),pointer :: io_domain
    integer,dimension(:),allocatable :: c
    integer,dimension(:),allocatable :: e
    if (is_unstructured) then
        allocate(c(size(fileobj%pelist)))
        allocate(e(size(fileobj%pelist)))
        io_domain => mpp_get_ug_io_domain(fileobj%domain)
        call mpp_get_ug_compute_domains(io_domain, &
                                        begin=c, &
                                        size=e)
        if (c(1) .ne. 1) then
            c(:) = c(:) - c(1) + 1
        endif
        call register_compressed_dimension(fileobj, &
                                           dim_name, &
                                           npes_corner=c, &
                                           npes_nelems=e)
        deallocate(c)
        deallocate(e)
    else
        call register_compressed_dimension(fileobj, &
                                           dim_name, &
                                           dim_size=dim_size)
    endif
end subroutine register_unstructured_dimension


end module fms_netcdf_unstructured_domain_io_mod
