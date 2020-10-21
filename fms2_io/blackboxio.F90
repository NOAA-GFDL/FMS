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

module blackboxio
use netcdf
use mpp_domains_mod
use fms_io_utils_mod
use netcdf_io_mod
use fms_netcdf_domain_io_mod
use fms_netcdf_unstructured_domain_io_mod
use mpp_mod, only: mpp_pe
use, intrinsic :: iso_fortran_env, only: error_unit, int32, int64, real32, real64
implicit none
private

integer, private :: fms2_ncchksz = -1 !< Chunksize (bytes) used in nc_open and nc_create

public :: blackboxio_init
public :: create_diskless_netcdf_file_wrap
public :: netcdf_save_restart_wrap2
public :: netcdf_restore_state_wrap
public :: create_diskless_domain_file
public :: save_domain_restart_wrap
public :: restore_domain_state_wrap
public :: create_diskless_unstructured_domain_file
public :: unstructured_write_restart_wrap


contains
!> @brief Accepts the namelist fms2_io_nml variables relevant to blackboxio 
subroutine blackboxio_init (chksz)
integer, intent(in) :: chksz
 fms2_ncchksz = chksz
end subroutine blackboxio_init


!> @brief Create a new file path.
!! @internal
subroutine get_new_filename(path, new_path, directory, timestamp, new_name)

  character(len=*), intent(in) :: path !< File path.
  character(len=*), intent(out) :: new_path !< New file path.
  character(len=*), intent(in), optional :: directory !< Directory
  character(len=*), intent(in), optional :: timestamp !< Time.
  character(len=*), intent(in), optional :: new_name !< New file basename.

  character(len=256) :: dir
  character(len=256) :: tstamp
  character(len=256) :: nname

  dir = ""
  if (present(directory)) then
    call string_copy(dir, trim(directory)//"/")
  endif
  tstamp = ""
  if (present(timestamp)) then
    call string_copy(tstamp, trim(timestamp)//".")
  endif
  call string_copy(nname, trim(path))
  if (present(new_name)) then
    call string_copy(nname, trim(new_name))
  endif
  call string_copy(new_path, trim(dir)//trim(tstamp)//trim(nname))
end subroutine get_new_filename


!> @brief Create a unique filename (poor man's version of mktemp).
!! @internal
subroutine tempfile(filename)

  character(len=*), intent(out) :: filename !< New unique filename.

  real :: numr
  integer :: numi

  do while(.true.)
    call random_number(numr)
    numi = transfer(numr, numi)
    numi = iand(numi, z'FFFFFF')
    write(filename, '(a,z6.6)') "tmp", numi
    if (.not. file_exists(filename)) then
      exit
    endif
  enddo
end subroutine tempfile


!> @brief Create a "diskless" netcdf file to act as a buffer to support our "register
!!        data to a file without knowing its name" legacy restart I/O workflow.
!! @return Flag telling whether the creation of the buffer was successful.
!! @internal
function create_diskless_netcdf_file(fileobj, pelist, path) &
  result(success)

  class(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  integer, dimension(:), intent(in), optional :: pelist !< List of ranks associated
                                                        !! with this file.  If not
                                                        !! provided, only the current
                                                        !! rank will be able to
                                                        !! act on the file.
  character(len=*), intent(in), optional :: path !< File path.
  logical :: success

  integer :: cmode
  integer :: err

  if (present(path)) then
    call string_copy(fileobj%path, path)
  else
    call tempfile(fileobj%path)
  endif
  fileobj%nc_format = "classic"
  fileobj%is_readonly = .false.
  if (present(pelist)) then
    allocate(fileobj%pelist(size(pelist)))
    fileobj%pelist(:) = pelist(:)
  else
    allocate(fileobj%pelist(1))
    fileobj%pelist(1) = mpp_pe()
  endif
  fileobj%io_root = fileobj%pelist(1)
  fileobj%is_root = mpp_pe() .eq. fileobj%io_root
  fileobj%is_restart = .true.
  fileobj%is_diskless = .true.
  cmode = ior(nf90_noclobber, nf90_classic_model)
  cmode = ior(cmode, nf90_diskless)
  if (fms2_ncchksz == -1) call error("create_diskless_netcdf_file :: fms2_ncchksz not set.")
  err = nf90_create(trim(fileobj%path), cmode, fileobj%ncid, chunksize=fms2_ncchksz)
  success = err .eq. nf90_noerr
  if (.not. success) then
    deallocate(fileobj%pelist)
    return
  endif
  allocate(fileobj%restart_vars(max_num_restart_vars))
  fileobj%num_restart_vars = 0
  allocate(fileobj%compressed_dims(max_num_compressed_dims))
  fileobj%num_compressed_dims = 0
end function create_diskless_netcdf_file


!> @brief Copy metadata from one file object to another.
!! @internal
subroutine copy_metadata(fileobj, new_fileobj)

  class(FmsNetcdfFile_t), intent(in), target :: fileobj !< File object.
  class(FmsNetcdfFile_t), intent(inout) :: new_fileobj !< New file object.

  integer :: err
  integer :: natt
  integer :: ndim
  integer :: varndim
  integer :: nvar
  character(len=nf90_max_name) :: n
  character(len=nf90_max_name) :: varname
  character(len=nf90_max_name), dimension(nf90_max_dims) :: dimnames
  integer, dimension(nf90_max_dims) :: dimlens
  integer :: xtype
  integer, dimension(nf90_max_var_dims) :: dimids
  integer, dimension(nf90_max_var_dims) :: d
  integer :: ulim_dimid
  integer :: varid
  integer :: i
  integer :: j
  integer :: k
  integer(kind=int32), dimension(:), allocatable :: buf_int
  real(kind=real32), dimension(:), allocatable :: buf_float
  real(kind=real64), dimension(:), allocatable :: buf_double

  if (fileobj%is_root .and. .not. new_fileobj%is_readonly) then
    !Copy global attributes to the new file.
    call set_netcdf_mode(fileobj%ncid, define_mode)
    call set_netcdf_mode(new_fileobj%ncid, define_mode)
    err = nf90_inquire(fileobj%ncid, nattributes=natt)
    call check_netcdf_code(err)
    do i = 1, natt
      err = nf90_inq_attname(fileobj%ncid, nf90_global, i, n)
      call check_netcdf_code(err)
      err = nf90_copy_att(fileobj%ncid, nf90_global, n, new_fileobj%ncid, nf90_global)
      call check_netcdf_code(err)
    enddo

    !Copy the dimensions to the new file.
    err = nf90_inquire(fileobj%ncid, ndimensions=ndim)
    call check_netcdf_code(err)
    err = nf90_inquire(fileobj%ncid, unlimiteddimid=ulim_dimid)
    call check_netcdf_code(err)
    do i = 1, ndim
      err = nf90_inquire_dimension(fileobj%ncid, i, dimnames(i), dimlens(i))
      call check_netcdf_code(err)
      if (i .eq. ulim_dimid) then
        err = nf90_def_dim(new_fileobj%ncid, dimnames(i), nf90_unlimited, dimids(i))
        ulim_dimid = dimids(i)
      else
        err = nf90_def_dim(new_fileobj%ncid, dimnames(i), dimlens(i), dimids(i))
      endif
      call check_netcdf_code(err)
    enddo

    !Copy the variables to the new file.
    err = nf90_inquire(fileobj%ncid, nvariables=nvar)
    call check_netcdf_code(err)
    do i = 1, nvar
      err = nf90_inquire_variable(fileobj%ncid, i, varname, xtype, varndim, d, natt)
      call check_netcdf_code(err)

      !Map to new dimension ids.
      do j = 1, varndim
        err = nf90_inquire_dimension(fileobj%ncid, d(j), n)
        call check_netcdf_code(err)
        do k = 1, ndim
          if (string_compare(n, dimnames(k))) then
            d(j) = dimids(k)
            exit
          endif
        enddo
      enddo

      !Define variable in new file.
      err = nf90_def_var(new_fileobj%ncid, varname, xtype, d(1:varndim), varid)
      call check_netcdf_code(err)

      !If the variable is an "axis", copy its data to the new file.
      if (varndim .eq. 1 .and. d(1) .ne. ulim_dimid) then
        do k = 1, ndim
          if (string_compare(varname, dimnames(k))) then
            call set_netcdf_mode(fileobj%ncid, data_mode)
            call set_netcdf_mode(new_fileobj%ncid, data_mode)
            if (xtype .eq. nf90_int) then
              allocate(buf_int(dimlens(k)))
              err = nf90_get_var(fileobj%ncid, i, buf_int)
              call check_netcdf_code(err)
              err = nf90_put_var(new_fileobj%ncid, varid, buf_int)
              deallocate(buf_int)
            elseif (xtype .eq. nf90_float) then
              allocate(buf_float(dimlens(k)))
              err = nf90_get_var(fileobj%ncid, i, buf_float)
              call check_netcdf_code(err)
              err = nf90_put_var(new_fileobj%ncid, varid, buf_float)
              deallocate(buf_float)
            elseif (xtype .eq. nf90_double) then
              allocate(buf_double(dimlens(k)))
              err = nf90_get_var(fileobj%ncid, i, buf_double)
              call check_netcdf_code(err)
              err = nf90_put_var(new_fileobj%ncid, varid, buf_double)
              deallocate(buf_double)
            else
              call error("this branch should not be reached.")
            endif
            call check_netcdf_code(err)
            call set_netcdf_mode(fileobj%ncid, define_mode)
            call set_netcdf_mode(new_fileobj%ncid, define_mode)
            exit
          endif
        enddo
      endif

      !Copy variable attributes to the new file.
      do j = 1, natt
        err = nf90_inq_attname(fileobj%ncid, i, j, n)
        call check_netcdf_code(err)
        err = nf90_copy_att(fileobj%ncid, i, n, new_fileobj%ncid, varid)
        call check_netcdf_code(err)
      enddo
    enddo
  endif

  if (new_fileobj%is_restart) then
    !Copy pointers to buffers (this is aliasing!).
    do i = 1, fileobj%num_restart_vars
      new_fileobj%restart_vars(i)%varname = fileobj%restart_vars(i)%varname
      if (associated(fileobj%restart_vars(i)%data0d)) then
        new_fileobj%restart_vars(i)%data0d => fileobj%restart_vars(i)%data0d
      elseif (associated(fileobj%restart_vars(i)%data1d)) then
        new_fileobj%restart_vars(i)%data1d => fileobj%restart_vars(i)%data1d
      elseif (associated(fileobj%restart_vars(i)%data2d)) then
        new_fileobj%restart_vars(i)%data2d => fileobj%restart_vars(i)%data2d
      elseif (associated(fileobj%restart_vars(i)%data3d)) then
        new_fileobj%restart_vars(i)%data3d => fileobj%restart_vars(i)%data3d
      elseif (associated(fileobj%restart_vars(i)%data4d)) then
        new_fileobj%restart_vars(i)%data4d => fileobj%restart_vars(i)%data4d
      else
        call error("this branch should not be reached.")
      endif
    enddo
    new_fileobj%num_restart_vars = fileobj%num_restart_vars
  endif

  !Copy compressed dimension metadata.
  do i = 1, fileobj%num_compressed_dims
    new_fileobj%compressed_dims(i)%dimname = fileobj%compressed_dims(i)%dimname
    k = size(fileobj%compressed_dims(i)%npes_corner)
    allocate(new_fileobj%compressed_dims(i)%npes_corner(k))
    allocate(new_fileobj%compressed_dims(i)%npes_nelems(k))
    do j = 1, k
      new_fileobj%compressed_dims(i)%npes_corner(j) = fileobj%compressed_dims(i)%npes_corner(j)
      new_fileobj%compressed_dims(i)%npes_nelems(j) = fileobj%compressed_dims(i)%npes_nelems(j)
    enddo
    new_fileobj%compressed_dims(i)%nelems = fileobj%compressed_dims(i)%nelems
  enddo
  new_fileobj%num_compressed_dims = fileobj%num_compressed_dims
end subroutine copy_metadata


!> @brief Make a copy of a file's metadata to support "intermediate restarts".
!! @internal
subroutine new_netcdf_file(fileobj, path, mode, new_fileobj, nc_format)

  class(FmsNetcdfFile_t), intent(in), target :: fileobj !< File object.
  character(len=*), intent(in) :: path !< Name of new file.
  character(len=*), intent(in) :: mode !< File mode.  Allowed values are:
                                       !! "read", "append", "write", or
                                       !! "overwrite".
  class(FmsNetcdfFile_t), intent(out) :: new_fileobj !< New file object.
  character(len=*), intent(in), optional :: nc_format !< Netcdf format that
                                                     !! new files are written
                                                     !! as.  Allowed values
                                                     !! are: "64bit", "classic",
                                                     !! or "netcdf4". Defaults to
                                                     !! "64bit".

  logical :: success

  !Open the new file.
  success = netcdf_file_open(new_fileobj, path, mode, nc_format, &
                             fileobj%pelist, fileobj%is_restart)
  if (.not. success) then
    call error("error opening file "//trim(path)//".")
  endif
  call copy_metadata(fileobj, new_fileobj)
end subroutine new_netcdf_file


!> @brief Wrapper to distinguish interfaces.
!! @return Flag telling whether the creation of the buffer was successful.
function create_diskless_netcdf_file_wrap(fileobj, pelist, path) &
  result(success)

  type(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  integer, dimension(:), intent(in), optional :: pelist !< List of ranks associated
                                                        !! with this file.  If not
                                                        !! provided, only the current
                                                        !! rank will be able to
                                                        !! act on the file.
  character(len=*), intent(in), optional :: path !< File path.
  logical :: success

  success = create_diskless_netcdf_file(fileobj, pelist, path)
end function create_diskless_netcdf_file_wrap


!> @brief Support for writing new restarts from a diskless file.
subroutine netcdf_save_restart_wrap2(fileobj, unlim_dim_level, directory, timestamp, &
                                     filename, nc_format)

  type(FmsNetcdfFile_t), intent(in), target :: fileobj !< File object.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  character(len=*), intent(in), optional :: directory !< Directory to write restart file to.
  character(len=*), intent(in), optional :: timestamp !< Model time.
  character(len=*), intent(in), optional :: filename !< New name for the file.
  character(len=*), intent(in), optional :: nc_format !< Netcdf format that
                                                     !! new files are written
                                                     !! as.  Allowed values
                                                     !! are: "64bit", "classic",
                                                     !! or "netcdf4". Defaults to
                                                     !! "64bit".

  character(len=256) :: new_name
  type(FmsNetcdfFile_t), target :: new_fileobj
  type(FmsNetcdfFile_t), pointer :: p
  logical :: close_new_file

  call get_new_filename(fileobj%path, new_name, directory, timestamp, filename)
  if (string_compare(fileobj%path, new_name)) then
    p => fileobj
    close_new_file = .false.
  else
    call new_netcdf_file(fileobj, new_name, "write", new_fileobj, nc_format)
    p => new_fileobj
    close_new_file = .true.
  endif
  call netcdf_save_restart(p, unlim_dim_level)
  if (close_new_file) then
    call netcdf_file_close(p)
  endif
end subroutine netcdf_save_restart_wrap2


!> @brief Loop through registered restart variables and read them from
!!        a netcdf file.
subroutine netcdf_restore_state_wrap(fileobj, unlim_dim_level, directory, timestamp, &
                                     filename)

  type(FmsNetcdfFile_t), intent(inout), target :: fileobj !< File object.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
  character(len=*), intent(in), optional :: directory !< Directory to write restart file to.
  character(len=*), intent(in), optional :: timestamp !< Model time.
  character(len=*), intent(in), optional :: filename !< New name for the file.

  character(len=256) :: new_name
  type(FmsNetcdfFile_t), target :: new_fileobj
  type(FmsNetcdfFile_t), pointer :: p
  logical :: close_new_file

  call get_new_filename(fileobj%path, new_name, directory, timestamp, filename)
  if (string_compare(fileobj%path, new_name)) then
    p => fileobj
    close_new_file = .false.
  else
    call new_netcdf_file(fileobj, new_name, "read", new_fileobj)
    p => new_fileobj
    close_new_file = .true.
  endif
  call netcdf_restore_state(p, unlim_dim_level)
  if (close_new_file) then
    call netcdf_file_close(p)
  endif
end subroutine netcdf_restore_state_wrap


!> @brief Create a "diskless" netcdf file to act as a buffer to support our "register
!!        data to a file without knowing its name" legacy restart I/O workflow.
!! @return Flag telling whether the creation of the buffer was successful.
function create_diskless_domain_file(fileobj, domain, path) &
  result(success)

  type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj !< File object.
  type(domain2d), intent(in) :: domain !< Two-dimensional domain.
  character(len=*), intent(in), optional :: path !< File path.
  logical :: success

  type(domain2d), pointer :: io_domain
  integer :: pelist_size
  integer, dimension(:), allocatable :: pelist

  io_domain => mpp_get_io_domain(domain)
  if (.not. associated(io_domain)) then
    call error("input domain does not have an io_domain.")
  endif
  pelist_size = mpp_get_domain_npes(io_domain)
  allocate(pelist(pelist_size))
  call mpp_get_pelist(io_domain, pelist)
  success = create_diskless_netcdf_file(fileobj, pelist, path)
  deallocate(pelist)
  fileobj%domain = domain
  allocate(fileobj%xdims(max_num_domain_decomposed_dims))
  fileobj%nx = 0
  allocate(fileobj%ydims(max_num_domain_decomposed_dims))
  fileobj%ny = 0
  call string_copy(fileobj%non_mangled_path, fileobj%path)
end function create_diskless_domain_file


!> @brief Make a copy of a file's metadata to support "intermediate restarts".
!! @internal
subroutine new_domain_file(fileobj, path, mode, new_fileobj, nc_format)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: path !< Name of new file.
  character(len=*), intent(in) :: mode !< File mode.  Allowed values are:
                                       !! "read", "append", "write", or "overwrite".
  type(FmsNetcdfDomainFile_t), intent(out) :: new_fileobj !< File object.
  character(len=*), intent(in), optional :: nc_format !< Netcdf format that
                                                     !! new files are written
                                                     !! as.  Allowed values
                                                     !! are: "64bit", "classic",
                                                     !! or "netcdf4". Defaults to
                                                     !! "64bit".

  logical :: success
  integer :: i

  success = open_domain_file(new_fileobj, path, mode, fileobj%domain, nc_format, &
                             fileobj%is_restart)
  if (.not. success) then
    call error("error opening file "//trim(path)//".")
  endif
  call copy_metadata(fileobj, new_fileobj)
  do i = 1, fileobj%nx
    call string_copy(new_fileobj%xdims(i)%varname, fileobj%xdims(i)%varname)
    new_fileobj%xdims(i)%pos = fileobj%xdims(i)%pos
  enddo
  new_fileobj%nx = fileobj%nx
  do i = 1, fileobj%ny
    call string_copy(new_fileobj%ydims(i)%varname, fileobj%ydims(i)%varname)
    new_fileobj%ydims(i)%pos = fileobj%ydims(i)%pos
  enddo
  new_fileobj%ny = fileobj%ny
end subroutine new_domain_file


!> @brief Loop through registered restart variables and write them to
!!        a netcdf file.
subroutine save_domain_restart_wrap(fileobj, unlim_dim_level, directory, timestamp, &
                                    filename, nc_format)

  type(FmsNetcdfDomainFile_t), intent(in), target :: fileobj !< File object.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.
  character(len=*), intent(in), optional :: directory !< Directory to write restart file to.
  character(len=*), intent(in), optional :: timestamp !< Model time.
  character(len=*), intent(in), optional :: filename !< New name for the file.
  character(len=*), intent(in), optional :: nc_format !< Netcdf format that
                                                     !! new files are written
                                                     !! as.  Allowed values
                                                     !! are: "64bit", "classic",
                                                     !! or "netcdf4". Defaults to
                                                     !! "64bit".

  character(len=256) :: new_name
  type(FmsNetcdfDomainFile_t), target :: new_fileobj
  type(FmsNetcdfDomainFile_t), pointer :: p
  logical :: close_new_file

  call get_new_filename(fileobj%non_mangled_path, new_name, directory, timestamp, filename)
  if (string_compare(fileobj%non_mangled_path, new_name)) then
    p => fileobj
    close_new_file = .false.
  else
    call new_domain_file(fileobj, new_name, "write", new_fileobj, nc_format)
    p => new_fileobj
    close_new_file = .true.
  endif
  call save_domain_restart(p, unlim_dim_level)
  if (close_new_file) then
    call close_domain_file(p)
  endif
end subroutine save_domain_restart_wrap


!> @brief Loop through registered restart variables and read them from
!!        a netcdf file.
subroutine restore_domain_state_wrap(fileobj, unlim_dim_level, directory, timestamp, &
                                     filename)

  type(FmsNetcdfDomainFile_t), intent(in), target :: fileobj !< File object.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.
  character(len=*), intent(in), optional :: directory !< Directory to write restart file to.
  character(len=*), intent(in), optional :: timestamp !< Model time.
  character(len=*), intent(in), optional :: filename !< New name for the file.

  character(len=256) :: new_name
  type(FmsNetcdfDomainFile_t), target :: new_fileobj
  type(FmsNetcdfDomainFile_t), pointer :: p
  logical :: close_new_file

  call get_new_filename(fileobj%non_mangled_path, new_name, directory, timestamp, filename)
  if (string_compare(fileobj%non_mangled_path, new_name)) then
    p => fileobj
    close_new_file = .false.
  else
    call new_domain_file(fileobj, new_name, "read", new_fileobj)
    p => new_fileobj
    close_new_file = .true.
  endif
  call restore_domain_state(p, unlim_dim_level)
  if (close_new_file) then
    call close_domain_file(p)
  endif
end subroutine restore_domain_state_wrap


!> @brief Create a "diskless" netcdf file to act as a buffer to support our "register
!!        data to a file without knowing its name" legacy restart I/O workflow.
!! @return Flag telling whether the creation of the buffer was successful.
function create_diskless_unstructured_domain_file(fileobj, domain, path) &
  result(success)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(inout) :: fileobj !< File object.
  type(domainug), intent(in) :: domain !< Two-dimensional domain.
  character(len=*), intent(in), optional :: path !< File path.
  logical :: success

  type(domainug), pointer :: io_domain
  integer :: pelist_size
  integer, dimension(:), allocatable :: pelist

  io_domain => mpp_get_ug_io_domain(domain)
  if (.not. associated(io_domain)) then
    call error("input domain does not have an io_domain.")
  endif
  pelist_size = mpp_get_ug_domain_npes(io_domain)
  allocate(pelist(pelist_size))
  call mpp_get_ug_domain_pelist(io_domain, pelist)
  success = create_diskless_netcdf_file(fileobj, pelist, path)
  deallocate(pelist)
  fileobj%domain = domain
  call string_copy(fileobj%non_mangled_path, fileobj%path)
end function create_diskless_unstructured_domain_file


!> @brief Make a copy of a file's metadata to support "intermediate restarts".
!! @internal
subroutine new_unstructured_domain_file(fileobj, path, mode, new_fileobj, nc_format)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: path !< Name of new file.
  character(len=*), intent(in) :: mode !< File mode.  Allowed values are:
                                       !! "read", "append", "write", or "overwrite."
  type(FmsNetcdfUnstructuredDomainFile_t), intent(out) :: new_fileobj !< New file object.
  character(len=*), intent(in), optional :: nc_format !< Netcdf format that
                                                     !! new files are written
                                                     !! as.  Allowed values
                                                     !! are: "64bit", "classic",
                                                     !! or "netcdf4". Defaults to
                                                     !! "64bit".

  logical :: success

  success = open_unstructured_domain_file(new_fileobj, path, mode, fileobj%domain, &
                                          nc_format, fileobj%is_restart)
  if (.not. success) then
    call error("error while opening file "//trim(path)//".")
  endif
  call copy_metadata(fileobj, new_fileobj)
end subroutine new_unstructured_domain_file


!> @brief Wrapper to distinguish interfaces.
subroutine unstructured_write_restart_wrap(fileobj, unlim_dim_level, directory, timestamp, &
                                           filename, nc_format)

  type(FmsNetcdfUnstructuredDomainFile_t), intent(in) :: fileobj !< File object.
  integer, intent(in), optional :: unlim_dim_level !< Unlimited dimension level.
  character(len=*), intent(in), optional :: directory !< Directory to write restart file to.
  character(len=*), intent(in), optional :: timestamp !< Model time.
  character(len=*), intent(in), optional :: filename !< New name for the file.
  character(len=*), intent(in), optional :: nc_format !< Netcdf format that
                                                     !! new files are written
                                                     !! as.  Allowed values
                                                     !! are: "64bit", "classic",
                                                     !! or "netcdf4". Defaults to
                                                     !! "64bit".

  character(len=256) :: new_name
  type(FmsNetcdfUnstructuredDomainFile_t) :: new_fileobj

  call get_new_filename(fileobj%non_mangled_path, new_name, directory, timestamp, filename)
  if (string_compare(fileobj%non_mangled_path, new_name)) then
    call netcdf_save_restart(fileobj, unlim_dim_level)
  else
    call new_unstructured_domain_file(fileobj, new_name, "write", new_fileobj, nc_format)
    call netcdf_save_restart(new_fileobj, unlim_dim_level)
    call close_unstructured_domain_file(new_fileobj)
  endif
end subroutine unstructured_write_restart_wrap


end module blackboxio
