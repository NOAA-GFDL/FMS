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
!> @file
!> @brief Routines for writing domain decomposed variables for the @ref write_data interface

!> @addtogroup fms2_io_mod
!> @{

!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_0d(fileobj, variable_name, vdata, unlim_dim_level, corner)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), intent(in) :: vdata !< Data that will
                                !! be written out
                                !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, intent(in), optional :: corner !< Array of starting
                                          !! indices describing
                                          !! where the data
                                          !! will be written to.

  call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner)

end subroutine domain_write_0d


!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_1d(fileobj, variable_name, vdata, unlim_dim_level, &
                           corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(in) :: vdata !< Data that will
                                              !! be written out
                                              !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(1), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(1), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)

end subroutine domain_write_1d


!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed" using mpp_gather.
subroutine domain_write_2d(fileobj, variable_name, vdata, unlim_dim_level, &
                           corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), contiguous, dimension(:,:), target, intent(in) :: vdata !< Data that will
                                                                    !! be written out
                                                                    !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(2), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(2), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer(kind=i4_kind), dimension(:,:), allocatable :: global_buf_i4_kind
  integer(kind=i8_kind), dimension(:,:), allocatable :: global_buf_i8_kind
  real(kind=r4_kind), dimension(:,:), allocatable :: global_buf_r4_kind
  real(kind=r8_kind), dimension(:,:), allocatable :: global_buf_r8_kind
  logical :: buffer_includes_halos
  type(domain2d), pointer :: io_domain
  integer :: isc
  integer :: isd
  integer :: jsc
  integer :: jsd
  integer :: xc_size
  integer :: xdim_index
  integer :: xpos
  integer :: ydim_index
  integer :: ypos
  integer :: yc_size
  integer(kind=i4_kind) :: fill_i4_kind !< Fill value of a i4_kind variable
  integer(kind=i8_kind) :: fill_i8_kind !< Fill value of a i8_kind variable
  real(kind=r4_kind) :: fill_r4_kind !< Fill value of a r4_kind variable
  real(kind=r8_kind) :: fill_r8_kind !< Fill value of a r8_kind variable
  class(*), dimension(:,:,:,:), pointer :: vdata_dummy    !< Vdata remapped as 4D
  integer :: xgmin !< Starting x index of the global io domain
  integer :: ygmin !< Ending y index of the global io domain
  integer :: xgsize, ygsize
  integer :: istart, jstart, iend, jend
  integer :: ioff, joff

  if (fileobj%use_netcdf_mpi) then
    call netcdf_mpi_write_2d(fileobj, variable_name, vdata, unlim_dim_level, corner, edge_lengths)
    return
  endif

  ! If the file is not domain decomposed, do not do anything in this routine
  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index, xpos, ypos)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)
    return
  endif

  if (xdim_index .ne. 1 .or. ydim_index .ne. 2) then
  ! This is a KLUDGE
  ! mpp_scatter assumes that the variable is (x,y), if that is not the case it remaps the data
  ! to a 4D array and calls domain_read_4d which does not use mpp_scatter yet
    vdata_dummy(1:size(vdata,1),1:size(vdata,2), 1:1, 1:1) => vdata(:,:)
    call domain_write_4d(fileobj, variable_name, vdata_dummy, unlim_dim_level)
    return
  endif

  ! Get the io domain from the fileobj
  io_domain => mpp_get_io_domain(fileobj%domain)

  ! Get some info about the domain:
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                        xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, buffer_includes_halos, &
                        msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))

  ! Get the global io domain:
  call mpp_get_global_domain(io_domain, xbegin=xgmin, xsize=xgsize, position=xpos)
  call mpp_get_global_domain(io_domain, ybegin=ygmin, ysize=ygsize, position=ypos)

  ! Root pe allocates room to gather the data and computes offsets for data placement
  if (fileobj%is_root) then
     ! Offsets used to place data in recv buffer
     ioff = 1-xgmin
     joff = 1-ygmin

     ! Allocate recv buffer for gather
     select type(vdata)
      type is (integer(kind=i4_kind))
        allocate(global_buf_i4_kind(xgsize, ygsize))
        global_buf_i4_kind = 0
        if (get_fill_value(fileobj, variable_name, fill_i4_kind, broadcast=.false.)) then
            global_buf_i4_kind = fill_i4_kind
        endif
      type is (integer(kind=i8_kind))
        allocate(global_buf_i8_kind(xgsize, ygsize))
        global_buf_i8_kind = 0
        if (get_fill_value(fileobj, variable_name, fill_i8_kind, broadcast=.false.)) then
            global_buf_i8_kind = fill_i8_kind
        endif
      type is (real(kind=r4_kind))
        allocate(global_buf_r4_kind(xgsize, ygsize))
        global_buf_r4_kind = 0.
        if (get_fill_value(fileobj, variable_name, fill_r4_kind, broadcast=.false.)) then
            global_buf_r4_kind = fill_r4_kind
        endif
      type is (real(kind=r8_kind))
        allocate(global_buf_r8_kind(xgsize, ygsize))
        global_buf_r8_kind = 0.
        if (get_fill_value(fileobj, variable_name, fill_r8_kind, broadcast=.false.)) then
            global_buf_r8_kind = fill_r8_kind
        endif
      class default
        call error("unsupported variable type: domain_write_2d_mpp_gather: file: " &
                 & //trim(fileobj%path)//" variable:"//trim(variable_name))
     end select
  endif

  ! Get the starting and indices of the compute domain relative to vdata (note that vdata start indices at 1 #Fortran)
  istart = 1
  jstart = 1

  ! If the buffer contains halos, get the portion of vdata with only the compute domain
  if (buffer_includes_halos) then
     istart = isc - isd + 1
     jstart = jsc - jsd + 1
  endif

  iend = istart + xc_size - 1
  jend = jstart + yc_size - 1

  ! Gather the data
  select type(vdata)
    type is (integer(kind=i4_kind))
      call mpp_gather(isc, isc+xc_size-1, jsc, jsc+yc_size-1, fileobj%pelist, vdata(istart:iend, jstart:jend), &
                      & global_buf_i4_kind, fileobj%is_root, ishift=ioff, jshift=joff)
    type is (integer(kind=i8_kind))
      call mpp_gather(isc, isc+xc_size-1, jsc, jsc+yc_size-1, fileobj%pelist, vdata(istart:iend, jstart:jend), &
                      & global_buf_i8_kind, fileobj%is_root, ishift=ioff, jshift=joff)
    type is (real(kind=r4_kind))
      call mpp_gather(isc, isc+xc_size-1, jsc, jsc+yc_size-1, fileobj%pelist, vdata(istart:iend, jstart:jend), &
                      & global_buf_r4_kind, fileobj%is_root, ishift=ioff, jshift=joff)
    type is (real(kind=r8_kind))
      call mpp_gather(isc, isc+xc_size-1, jsc, jsc+yc_size-1, fileobj%pelist, vdata(istart:iend, jstart:jend), &
                      & global_buf_r8_kind, fileobj%is_root, ishift=ioff, jshift=joff)
    class default
      call error("unsupported variable type: domain_write_2d_mpp_gather: file: " &
               & //trim(fileobj%path)//" variable:"// trim(variable_name))
  end select

  ! Root pe writes out the data
  if (fileobj%is_root) then
    select type(vdata)
      type is (integer(kind=i4_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_i4_kind, &
                               & unlim_dim_level=unlim_dim_level)
      type is (integer(kind=i8_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_i8_kind, &
                               & unlim_dim_level=unlim_dim_level)
      type is (real(kind=r4_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_r4_kind, &
                               & unlim_dim_level=unlim_dim_level)
      type is (real(kind=r8_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_r8_kind, &
                               & unlim_dim_level=unlim_dim_level)
      class default
        call error("unsupported variable type: domain_write_2d_mpp_gather: file: " &
                 & //trim(fileobj%path)//" variable:"// trim(variable_name))
    end select
  endif

end subroutine domain_write_2d


!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed" using mpp_gather.
subroutine domain_write_3d(fileobj, variable_name, vdata, unlim_dim_level, &
                                      corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), contiguous, dimension(:,:,:), target, intent(in) :: vdata !< Data that will
                                                                      !! be written out
                                                                      !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(3), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(3), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer(kind=i4_kind), dimension(:,:,:), allocatable :: global_buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:), allocatable :: global_buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:), allocatable :: global_buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:), allocatable :: global_buf_r8_kind
  logical :: buffer_includes_halos
  type(domain2d), pointer :: io_domain
  integer :: isc
  integer :: isd
  integer :: jsc
  integer :: jsd
  integer :: xc_size
  integer :: xdim_index
  integer :: xpos
  integer :: ydim_index
  integer :: ypos
  integer :: yc_size
  integer(kind=i4_kind) :: fill_i4_kind !< Fill value of a i4_kind variable
  integer(kind=i8_kind) :: fill_i8_kind !< Fill value of a i8_kind variable
  real(kind=r4_kind) :: fill_r4_kind !< Fill value of a r4_kind variable
  real(kind=r8_kind) :: fill_r8_kind !< Fill value of a r8_kind variable
  class(*), dimension(:,:,:,:), pointer :: vdata_dummy    !< Vdata remapped as 4D
  integer :: xgmin !< Starting x index of the global io domain
  integer :: ygmin !< Ending y index of the global io domain
  integer :: xgsize, ygsize
  integer :: istart, jstart, iend, jend
  integer :: ioff, joff

  if (fileobj%use_netcdf_mpi) then
    call netcdf_mpi_write_3d(fileobj, variable_name, vdata, unlim_dim_level, corner, edge_lengths)
    return
  endif

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                            xdim_index, ydim_index, xpos, ypos)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)
    return
  endif

  if (xdim_index .ne. 1 .or. ydim_index .ne. 2) then
  ! This is a KLUDGE
  ! mpp_scatter assumes that the variable is (x,y), if that is not the case it remaps the data
  ! to a 4D array and calls domain_read_4d which does not use mpp_scatter yet
    vdata_dummy(1:size(vdata,1),1:size(vdata,2), 1:size(vdata,3), 1:1) => vdata(:,:,:)
    call domain_write_4d(fileobj, variable_name, vdata_dummy, unlim_dim_level)
    return
  endif

  ! Get the io domain from the fileobj
  io_domain => mpp_get_io_domain(fileobj%domain)

  ! Get some info about the domain:
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, buffer_includes_halos, &
                      msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))

  ! Get the global io domain:
  call mpp_get_global_domain(io_domain, xbegin=xgmin, xsize=xgsize, position=xpos)
  call mpp_get_global_domain(io_domain, ybegin=ygmin, ysize=ygsize, position=ypos)

  ! Root pe allocates room to gather the data and computes offsets for data placement
  if (fileobj%is_root) then
     ! Offsets used to place data in recv buffer
     ioff = 1-xgmin
     joff = 1-ygmin

     ! Allocate recv buffer for gather
     select type(vdata)
      type is (integer(kind=i4_kind))
        allocate(global_buf_i4_kind(xgsize, ygsize, size(vdata,3)))
        global_buf_i4_kind = 0
        if (get_fill_value(fileobj, variable_name, fill_i4_kind, broadcast=.false.)) then
            global_buf_i4_kind = fill_i4_kind
        endif
      type is (integer(kind=i8_kind))
        allocate(global_buf_i8_kind(xgsize, ygsize, size(vdata,3)))
        global_buf_i8_kind = 0
        if (get_fill_value(fileobj, variable_name, fill_i8_kind, broadcast=.false.)) then
            global_buf_i8_kind = fill_i8_kind
        endif
      type is (real(kind=r4_kind))
        allocate(global_buf_r4_kind(xgsize, ygsize, size(vdata,3)))
        global_buf_r4_kind = 0.
        if (get_fill_value(fileobj, variable_name, fill_r4_kind, broadcast=.false.)) then
            global_buf_r4_kind = fill_r4_kind
        endif
      type is (real(kind=r8_kind))
        allocate(global_buf_r8_kind(xgsize, ygsize, size(vdata,3)))
        global_buf_r8_kind = 0.
        if (get_fill_value(fileobj, variable_name, fill_r8_kind, broadcast=.false.)) then
            global_buf_r8_kind = fill_r8_kind
        endif
      class default
        call error("unsupported variable type: domain_write_3d_mpp_gather: file: " &
                  & //trim(fileobj%path)//" variable:"//trim(variable_name))
     end select
  endif

  ! Get the starting and indices of the compute domain relative to vdata(note that vdata start indices at 1 #Fortran)
  istart = 1
  jstart = 1

  ! If the buffer contains halos, get the portion of vdata with only the compute domain
  if (buffer_includes_halos) then
     istart = isc - isd + 1
     jstart = jsc - jsd + 1
  endif

  iend = istart + xc_size - 1
  jend = jstart + yc_size - 1

  ! Get offsets for buffer
  ioff = 1-xgmin
  joff = 1-ygmin

  ! Gather the data
  select type(vdata)
    type is (integer(kind=i4_kind))
      call mpp_gather(isc, isc+xc_size-1, jsc, jsc+yc_size-1, size(vdata,3), fileobj%pelist, &
                  & vdata(istart:iend, jstart:jend,:), global_buf_i4_kind, fileobj%is_root, ishift=ioff, jshift=joff)
    type is (integer(kind=i8_kind))
      call mpp_gather(isc, isc+xc_size-1, jsc, jsc+yc_size-1, size(vdata,3), fileobj%pelist, &
                  & vdata(istart:iend, jstart:jend,:), global_buf_i8_kind, fileobj%is_root, ishift=ioff, jshift=joff)
    type is (real(kind=r4_kind))
      call mpp_gather(isc, isc+xc_size-1, jsc, jsc+yc_size-1, size(vdata,3), fileobj%pelist, &
                  & vdata(istart:iend, jstart:jend,:), global_buf_r4_kind, fileobj%is_root, ishift=ioff, jshift=joff)
    type is (real(kind=r8_kind))
      call mpp_gather(isc, isc+xc_size-1, jsc, jsc+yc_size-1, size(vdata,3), fileobj%pelist, &
                  & vdata(istart:iend, jstart:jend,:), global_buf_r8_kind, fileobj%is_root, ishift=ioff, jshift=joff)
    class default
      call error("unsupported variable type: domain_write_3d_mpp_gather: file: " &
                  & //trim(fileobj%path)//" variable:"//trim(variable_name))
  end select

  ! Root pe writes out the data
  if (fileobj%is_root) then
    select type(vdata)
      type is (integer(kind=i4_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_i4_kind, &
                               unlim_dim_level=unlim_dim_level)
      type is (integer(kind=i8_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_i8_kind, &
                               unlim_dim_level=unlim_dim_level)
      type is (real(kind=r4_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_r4_kind, &
                               unlim_dim_level=unlim_dim_level)
      type is (real(kind=r8_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_r8_kind, &
                               unlim_dim_level=unlim_dim_level)
      class default
        call error("unsupported variable type: domain_write_3d_mpp_gather: file: " &
                 & //trim(fileobj%path)//" variable:"//trim(variable_name))
    end select
  endif

end subroutine domain_write_3d


!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_4d(fileobj, variable_name, vdata, unlim_dim_level, &
                           corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(in) :: vdata !< Data that will
                                                    !! be written out
                                                    !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(4), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(4), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer(kind=i4_kind), dimension(:,:,:,:), allocatable :: buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:,:), allocatable :: buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:,:), allocatable :: buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:,:), allocatable :: buf_r8_kind
  logical :: buffer_includes_halos
  integer, dimension(4) :: c
  integer, dimension(4) :: e
  integer(kind=i4_kind), dimension(:,:,:,:), allocatable :: global_buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:,:), allocatable :: global_buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:,:), allocatable :: global_buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:,:), allocatable :: global_buf_r8_kind
  integer :: i
  type(domain2d), pointer :: io_domain
  integer :: isc
  integer :: isd
  integer :: jsc
  integer :: jsd
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_iec
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_jcsize
  integer, dimension(:), allocatable :: pe_jec
  integer, dimension(:), allocatable :: pe_jsc
  integer :: xc_size
  integer :: xdim_index
  integer :: xpos
  integer :: ydim_index
  integer :: ypos
  integer :: yc_size
  integer(kind=i4_kind) :: fill_i4_kind !< Fill value of a i4_kind variable
  integer(kind=i8_kind) :: fill_i8_kind !< Fill value of a i8_kind variable
  real(kind=r4_kind) :: fill_r4_kind !< Fill value of a r4_kind variable
  real(kind=r8_kind) :: fill_r8_kind !< Fill value of a r8_kind variable
  integer :: xgmax !< Ending x index of the global io domain
  integer :: xgmin !< Starting x index of the global io domain
  integer :: ygmax !< Ending y index of the global io domain
  integer :: ygmin !< Ending y index of the global io domain

  if (fileobj%use_netcdf_mpi) then
    call netcdf_mpi_write_4d(fileobj, variable_name, vdata, unlim_dim_level, corner, edge_lengths)
    return
  endif

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index, xpos, ypos)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, &
                      buffer_includes_halos, msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root gathers the data and writes it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_iec(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xend=pe_iec, xsize=pe_icsize, &
                                 position=xpos)
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jec(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, ybegin=pe_jsc, yend=pe_jec, ysize=pe_jcsize, &
                                 position=ypos)

    !< Determine the size of the global io domain
    call mpp_get_global_domain(io_domain, ybegin=ygmin, yend=ygmax, position=ypos)
    call mpp_get_global_domain(io_domain, xbegin=xgmin, xend=xgmax, position=xpos)

    !Write the out the data.
    !< Set e to equal the size of the global io domain
    e(xdim_index) = xgmax - xgmin + 1
    e(ydim_index) = ygmax - ygmin + 1

    !< Allocate a global buffer, get the fill value if it exists in the file, and initialize
    !! the buffer to the fill value
    select type(vdata)
      type is (integer(kind=i4_kind))
        call allocate_array(global_buf_i4_kind, e)
        global_buf_i4_kind = 0
        if (get_fill_value(fileobj, variable_name, fill_i4_kind, broadcast=.false.)) then
            global_buf_i4_kind = fill_i4_kind
        endif
      type is (integer(kind=i8_kind))
        call allocate_array(global_buf_i8_kind, e)
        global_buf_i8_kind = 0
        if (get_fill_value(fileobj, variable_name, fill_i8_kind, broadcast=.false.)) then
            global_buf_i8_kind = fill_i8_kind
        endif
      type is (real(kind=r4_kind))
        call allocate_array(global_buf_r4_kind, e)
        global_buf_r4_kind = 0.
        if (get_fill_value(fileobj, variable_name, fill_r4_kind, broadcast=.false.)) then
            global_buf_r4_kind = fill_r4_kind
        endif
      type is (real(kind=r8_kind))
        call allocate_array(global_buf_r8_kind, e)
        global_buf_r8_kind = 0.
        if (get_fill_value(fileobj, variable_name, fill_r8_kind, broadcast=.false.)) then
            global_buf_r8_kind = fill_r8_kind
        endif
      class default
        call error("unsupported variable type: domain_write_4d: file: "//trim(fileobj%path)//" variable:"// &
                 & trim(variable_name))
    end select

    do i = 1, size(fileobj%pelist)
      !< Set c relative to the starting global io domain index
      c(xdim_index) = pe_isc(i) - xgmin + 1
      c(ydim_index) = pe_jsc(i) - ygmin + 1
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=i4_kind))
          call allocate_array(buf_i4_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call get_array_section(buf_i4_kind, vdata, c, e)
            c(xdim_index) = pe_isc(i) - xgmin + 1
            c(ydim_index) = pe_jsc(i) - ygmin + 1
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_i4_kind, size(buf_i4_kind), fileobj%pelist(i), block=.true.)
          endif
          !Put local data into the global buffer.
          call put_array_section(buf_i4_kind, global_buf_i4_kind, c, e)
          deallocate(buf_i4_kind)
        type is (integer(kind=i8_kind))
          call allocate_array(buf_i8_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call get_array_section(buf_i8_kind, vdata, c, e)
            c(xdim_index) = pe_isc(i) - xgmin + 1
            c(ydim_index) = pe_jsc(i) - ygmin + 1
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_i8_kind, size(buf_i8_kind), fileobj%pelist(i), block=.true.)
          endif
          !Put local data into the global buffer.
          call put_array_section(buf_i8_kind, global_buf_i8_kind, c, e)
          deallocate(buf_i8_kind)
        type is (real(kind=r4_kind))
          call allocate_array(buf_r4_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call get_array_section(buf_r4_kind, vdata, c, e)
            c(xdim_index) = pe_isc(i) - xgmin + 1
            c(ydim_index) = pe_jsc(i) - ygmin + 1
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_r4_kind, size(buf_r4_kind), fileobj%pelist(i), block=.true.)
          endif
          !Put local data into the global buffer.
          call put_array_section(buf_r4_kind, global_buf_r4_kind, c, e)
          deallocate(buf_r4_kind)
        type is (real(kind=r8_kind))
          call allocate_array(buf_r8_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call get_array_section(buf_r8_kind, vdata, c, e)
            c(xdim_index) = pe_isc(i) - xgmin + 1
            c(ydim_index) = pe_jsc(i) - ygmin + 1
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_r8_kind, size(buf_r8_kind), fileobj%pelist(i), block=.true.)
          endif
          !Put local data into the global buffer.
          call put_array_section(buf_r8_kind, global_buf_r8_kind, c, e)
          deallocate(buf_r8_kind)
      end select
    enddo
    deallocate(pe_isc)
    deallocate(pe_iec)
    deallocate(pe_icsize)
    deallocate(pe_jsc)
    deallocate(pe_jec)
    deallocate(pe_jcsize)

    !Write the out the data.
    select type(vdata)
      type is (integer(kind=i4_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_i4_kind, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(global_buf_i4_kind)
      type is (integer(kind=i8_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_i8_kind, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(global_buf_i8_kind)
      type is (real(kind=r4_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_r4_kind, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(global_buf_r4_kind)
      type is (real(kind=r8_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_r8_kind, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(global_buf_r8_kind)
    end select
  else
    if (buffer_includes_halos) then
      c(xdim_index) = isc - isd + 1
      c(ydim_index) = jsc - jsd + 1
    endif
    e(xdim_index) = xc_size
    e(ydim_index) = yc_size
    select type(vdata)
      type is (integer(kind=i4_kind))
        call allocate_array(buf_i4_kind, e)
        call get_array_section(buf_i4_kind, vdata, c, e)
        call mpp_send(buf_i4_kind, size(buf_i4_kind), fileobj%io_root)
        call mpp_sync_self(check=event_send)
        deallocate(buf_i4_kind)
      type is (integer(kind=i8_kind))
        call allocate_array(buf_i8_kind, e)
        call get_array_section(buf_i8_kind, vdata, c, e)
        call mpp_send(buf_i8_kind, size(buf_i8_kind), fileobj%io_root)
        call mpp_sync_self(check=event_send)
        deallocate(buf_i8_kind)
      type is (real(kind=r4_kind))
        call allocate_array(buf_r4_kind, e)
        call get_array_section(buf_r4_kind, vdata, c, e)
        call mpp_send(buf_r4_kind, size(buf_r4_kind), fileobj%io_root)
        call mpp_sync_self(check=event_send)
        deallocate(buf_r4_kind)
      type is (real(kind=r8_kind))
        call allocate_array(buf_r8_kind, e)
        call get_array_section(buf_r8_kind, vdata, c, e)
        call mpp_send(buf_r8_kind, size(buf_r8_kind), fileobj%io_root)
        call mpp_sync_self(check=event_send)
        deallocate(buf_r8_kind)
      class default
        call error("unsupported variable type: domain_write_4d: file: "//trim(fileobj%path)//" variable:"// &
                 & trim(variable_name))
    end select
  endif
end subroutine domain_write_4d


!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_5d(fileobj, variable_name, vdata, unlim_dim_level, &
                           corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(in) :: vdata !< Data that will
                                                      !! be written out
                                                      !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, dimension(5), intent(in), optional :: corner !< Array of starting
                                                        !! indices describing
                                                        !! where the data
                                                        !! will be written to.
  integer, dimension(5), intent(in), optional :: edge_lengths !< The number of
                                                              !! elements that
                                                              !! will be written
                                                              !! in each dimension.

  integer(kind=i4_kind), dimension(:,:,:,:,:), allocatable :: buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:,:,:), allocatable :: buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:,:,:), allocatable :: buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:,:,:), allocatable :: buf_r8_kind
  logical :: buffer_includes_halos
  integer, dimension(5) :: c
  integer, dimension(5) :: e
  integer(kind=i4_kind), dimension(:,:,:,:,:), allocatable :: global_buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:,:,:), allocatable :: global_buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:,:,:), allocatable :: global_buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:,:,:), allocatable :: global_buf_r8_kind
  integer :: i
  type(domain2d), pointer :: io_domain
  integer :: isc
  integer :: isd
  integer :: jsc
  integer :: jsd
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_iec
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_jcsize
  integer, dimension(:), allocatable :: pe_jec
  integer, dimension(:), allocatable :: pe_jsc
  integer :: xc_size
  integer :: xdim_index
  integer :: xpos
  integer :: ydim_index
  integer :: ypos
  integer :: yc_size
  integer(kind=i4_kind) :: fill_i4_kind !< Fill value of a i4_kind variable
  integer(kind=i8_kind) :: fill_i8_kind !< Fill value of a i8_kind variable
  real(kind=r4_kind) :: fill_r4_kind !< Fill value of a r4_kind variable
  real(kind=r8_kind) :: fill_r8_kind !< Fill value of a r8_kind variable
  integer :: xgmax !< Ending x index of the global io domain
  integer :: xgmin !< Starting x index of the global io domain
  integer :: ygmax !< Ending y index of the global io domain
  integer :: ygmin !< Ending y index of the global io domain

  if (fileobj%use_netcdf_mpi) then
    call netcdf_mpi_write_5d(fileobj, variable_name, vdata, unlim_dim_level, corner, edge_lengths)
    return
  endif

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index, xpos, ypos)) then
    call compressed_write(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, &
                      buffer_includes_halos, msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root gathers the data and writes it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_iec(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xend=pe_iec, xsize=pe_icsize, &
                                 position=xpos)
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jec(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, ybegin=pe_jsc, yend=pe_jec, ysize=pe_jcsize, &
                                 position=ypos)

    !< Determine the size of the global io domain
    call mpp_get_global_domain(io_domain, ybegin=ygmin, yend=ygmax, position=ypos)
    call mpp_get_global_domain(io_domain, xbegin=xgmin, xend=xgmax, position=xpos)

    !Write the out the data.
    !< Set e to equal the size of the global io domain
    e(xdim_index) = xgmax - xgmin + 1
    e(ydim_index) = ygmax - ygmin + 1

    !< Allocate a global buffer, get the fill value if it exists in the file, and initialize
    !! the buffer to the fill value
    select type(vdata)
      type is (integer(kind=i4_kind))
        call allocate_array(global_buf_i4_kind, e)
        global_buf_i4_kind = 0
        if (get_fill_value(fileobj, variable_name, fill_i4_kind, broadcast=.false.)) then
            global_buf_i4_kind = fill_i4_kind
        endif
      type is (integer(kind=i8_kind))
        call allocate_array(global_buf_i8_kind, e)
        global_buf_i8_kind = 0
        if (get_fill_value(fileobj, variable_name, fill_i8_kind, broadcast=.false.)) then
            global_buf_i8_kind = fill_i8_kind
        endif
      type is (real(kind=r4_kind))
        call allocate_array(global_buf_r4_kind, e)
        global_buf_r4_kind = 0.
        if (get_fill_value(fileobj, variable_name, fill_r4_kind, broadcast=.false.)) then
            global_buf_r4_kind = fill_r4_kind
        endif
      type is (real(kind=r8_kind))
        call allocate_array(global_buf_r8_kind, e)
        global_buf_r8_kind = 0.
        if (get_fill_value(fileobj, variable_name, fill_r8_kind, broadcast=.false.)) then
            global_buf_r8_kind = fill_r8_kind
        endif
      class default
        call error("unsupported variable type: domain_write_5d: file: "//trim(fileobj%path)//" variable:"// &
                 & trim(variable_name))
    end select

    do i = 1, size(fileobj%pelist)
      !< Set c relative to the starting global io domain index
      c(xdim_index) = pe_isc(i) - xgmin + 1
      c(ydim_index) = pe_jsc(i) - ygmin + 1
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=i4_kind))
          call allocate_array(buf_i4_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call get_array_section(buf_i4_kind, vdata, c, e)
            c(xdim_index) = pe_isc(i) - xgmin + 1
            c(ydim_index) = pe_jsc(i) - ygmin + 1
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_i4_kind, size(buf_i4_kind), fileobj%pelist(i), block=.true.)
          endif
          !Put local data into the global buffer.
          call put_array_section(buf_i4_kind, global_buf_i4_kind, c, e)
          deallocate(buf_i4_kind)
        type is (integer(kind=i8_kind))
          call allocate_array(buf_i8_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call get_array_section(buf_i8_kind, vdata, c, e)
            c(xdim_index) = pe_isc(i) - xgmin + 1
            c(ydim_index) = pe_jsc(i) - ygmin + 1
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_i8_kind, size(buf_i8_kind), fileobj%pelist(i), block=.true.)
          endif
          !Put local data into the global buffer.
          call put_array_section(buf_i8_kind, global_buf_i8_kind, c, e)
          deallocate(buf_i8_kind)
        type is (real(kind=r4_kind))
          call allocate_array(buf_r4_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call get_array_section(buf_r4_kind, vdata, c, e)
            c(xdim_index) = pe_isc(i) - xgmin + 1
            c(ydim_index) = pe_jsc(i) - ygmin + 1
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_r4_kind, size(buf_r4_kind), fileobj%pelist(i), block=.true.)
          endif
          !Put local data into the global buffer.
          call put_array_section(buf_r4_kind, global_buf_r4_kind, c, e)
          deallocate(buf_r4_kind)
        type is (real(kind=r8_kind))
          call allocate_array(buf_r8_kind, e)
          !Get the data for fileobj%pelist(i)'s portion of the compute domain.
          if (i .eq. 1) then
            !Root rank gets the data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call get_array_section(buf_r8_kind, vdata, c, e)
            c(xdim_index) = pe_isc(i) - xgmin + 1
            c(ydim_index) = pe_jsc(i) - ygmin + 1
          else
            !Receive data from non-root ranks.
            call mpp_recv(buf_r8_kind, size(buf_r8_kind), fileobj%pelist(i), block=.true.)
          endif
          !Put local data into the global buffer.
          call put_array_section(buf_r8_kind, global_buf_r8_kind, c, e)
          deallocate(buf_r8_kind)
      end select
    enddo
    deallocate(pe_isc)
    deallocate(pe_iec)
    deallocate(pe_icsize)
    deallocate(pe_jsc)
    deallocate(pe_jec)
    deallocate(pe_jcsize)

    !Write the out the data.
    select type(vdata)
      type is (integer(kind=i4_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_i4_kind, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(global_buf_i4_kind)
      type is (integer(kind=i8_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_i8_kind, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(global_buf_i8_kind)
      type is (real(kind=r4_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_r4_kind, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(global_buf_r4_kind)
      type is (real(kind=r8_kind))
        call netcdf_write_data(fileobj, variable_name, global_buf_r8_kind, &
                               unlim_dim_level=unlim_dim_level)
        deallocate(global_buf_r8_kind)
    end select
  else
    if (buffer_includes_halos) then
      c(xdim_index) = isc - isd + 1
      c(ydim_index) = jsc - jsd + 1
    endif
    e(xdim_index) = xc_size
    e(ydim_index) = yc_size
    select type(vdata)
      type is (integer(kind=i4_kind))
        call allocate_array(buf_i4_kind, e)
        call get_array_section(buf_i4_kind, vdata, c, e)
        call mpp_send(buf_i4_kind, size(buf_i4_kind), fileobj%io_root)
        call mpp_sync_self(check=event_send)
        deallocate(buf_i4_kind)
      type is (integer(kind=i8_kind))
        call allocate_array(buf_i8_kind, e)
        call get_array_section(buf_i8_kind, vdata, c, e)
        call mpp_send(buf_i8_kind, size(buf_i8_kind), fileobj%io_root)
        call mpp_sync_self(check=event_send)
        deallocate(buf_i8_kind)
      type is (real(kind=r4_kind))
        call allocate_array(buf_r4_kind, e)
        call get_array_section(buf_r4_kind, vdata, c, e)
        call mpp_send(buf_r4_kind, size(buf_r4_kind), fileobj%io_root)
        call mpp_sync_self(check=event_send)
        deallocate(buf_r4_kind)
      type is (real(kind=r8_kind))
        call allocate_array(buf_r8_kind, e)
        call get_array_section(buf_r8_kind, vdata, c, e)
        call mpp_send(buf_r8_kind, size(buf_r8_kind), fileobj%io_root)
        call mpp_sync_self(check=event_send)
        deallocate(buf_r8_kind)
      class default
        call error("unsupported variable type: domain_write_5d: file: "//trim(fileobj%path)//" variable:"// &
                 & trim(variable_name))
    end select
  endif
end subroutine domain_write_5d

!> @brief Write 2D data using NetCDF MPI
subroutine netcdf_mpi_write_2d(fileobj, variable_name, vdata, unlim_dim_level, corner, edge_lengths)
  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in) :: vdata !< Data that will be written out to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited dimension.
  integer, dimension(2), intent(in), optional :: corner !< Array of starting indices describing where the data
                                                        !! will be written to.
  integer, dimension(2), intent(in), optional :: edge_lengths !< The number of elements that will be written
                                                              !! in each dimension.

  integer :: xdim_index
  integer :: ydim_index
  integer :: xpos
  integer :: ypos

  integer :: isd !< Starting "x" index of the data domain
  integer :: isc !< Starting "x" index of the compute domain
  integer :: jsd !< Starting "y" index of the data domain
  integer :: jsc !< Starting "y" index of the compute domain
  integer :: xc_size !< Size of the compute domain
  integer :: yc_size !< Size of the compute domain
  logical :: buffer_includes_halos !< True if "vdata" is the size of the data domain
                                   !! (i.e it includes halos (which are not written))

  integer, dimension(3) :: c !< Indices of the corners for each dimension
  integer, dimension(3) :: e !< Size of the data for each dimension
  integer :: varid
  integer :: unlim_dim_index

  integer :: i0, i1
  integer :: j0, j1
  integer :: err

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., xdim_index, ydim_index, xpos, ypos)) then
    call compressed_write(fileobj, variable_name, vdata, unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)
    return
  endif

  !! Get some more info about the variable
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, &
                      buffer_includes_halos, msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))

  c = 1
  e = [shape(vdata), 1]

  c(xdim_index) = isc
  c(ydim_index) = jsc
  e(xdim_index) = xc_size
  e(ydim_index) = yc_size

  if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, broadcast=.false.)
      if (unlim_dim_index .ne. 3) then
        call error("unlimited dimension must be the slowest varying dimension in variable: "//trim(variable_name))
      endif
      c(unlim_dim_index) = unlim_dim_level
  endif

  i0 = 1
  j0 = 1
  select case (xdim_index)
    case (1)
      ! X is the first dimension
      if (buffer_includes_halos) then
        i0 = isc - isd + 1
        j0 = jsc - jsd + 1
      endif
      i1 = i0 + xc_size - 1
      j1 = j0 + yc_size - 1
    case (2)
      ! X is the second dimension
      if (buffer_includes_halos) then
        i0 = jsc - jsd + 1
        j0 = isc - isd + 1
      endif
      i1 = i0 + yc_size - 1
      j1 = j0 + xc_size - 1
  end select

  varid = get_variable_id(fileobj%ncid, trim(variable_name), &
    msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))

  select type(vdata)
    type is (integer(kind=i4_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1), &
        start=c, count=e)
    type is (integer(kind=i8_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1), &
        start=c, count=e)
    type is (real(kind=r4_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1), &
        start=c, count=e)
    type is (real(kind=r8_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1), &
        start=c, count=e)
  end select

  call check_netcdf_code(err, "Failed to write variable: "//trim(variable_name))
end subroutine netcdf_mpi_write_2d

!> @brief Write 3D data using NetCDF MPI
subroutine netcdf_mpi_write_3d(fileobj, variable_name, vdata, unlim_dim_level, corner, edge_lengths)
  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in) :: vdata !< Data that will be written out to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited dimension.
  integer, dimension(3), intent(in), optional :: corner !< Array of starting indices describing where the data
                                                        !! will be written to.
  integer, dimension(3), intent(in), optional :: edge_lengths !< The number of elements that will be written
                                                              !! in each dimension.

  integer :: xdim_index
  integer :: ydim_index
  integer :: xpos
  integer :: ypos

  integer :: isd !< Starting "x" index of the data domain
  integer :: isc !< Starting "x" index of the compute domain
  integer :: jsd !< Starting "y" index of the data domain
  integer :: jsc !< Starting "y" index of the compute domain
  integer :: xc_size !< Size of the compute domain
  integer :: yc_size !< Size of the compute domain
  logical :: buffer_includes_halos !< True if "vdata" is the size of the data domain
                                   !! (i.e it includes halos (which are not written))

  integer, dimension(4) :: c !< Indices of the corners for each dimension
  integer, dimension(4) :: e !< Size of the data for each dimension
  integer :: varid
  integer :: unlim_dim_index

  integer :: i0, i1
  integer :: j0, j1
  integer :: k0, k1
  integer :: err

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., xdim_index, ydim_index, xpos, ypos)) then
    call compressed_write(fileobj, variable_name, vdata, unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)
    return
  endif

  !! Get some more info about the variable
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, &
                      buffer_includes_halos, msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))

  c = 1
  e = [shape(vdata), 1]

  c(xdim_index) = isc
  c(ydim_index) = jsc
  e(xdim_index) = xc_size
  e(ydim_index) = yc_size

  if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, broadcast=.false.)
      if (unlim_dim_index .ne. 4) then
        call error("unlimited dimension must be the slowest varying dimension in variable: "//trim(variable_name))
      endif
      c(unlim_dim_index) = unlim_dim_level
  endif

  i0 = 1
  j0 = 1
  k0 = 1

  i1 = size(vdata, 1)
  j1 = size(vdata, 2)
  k1 = size(vdata, 3)

  select case (xdim_index)
    case (1)
      ! X is the first dimension
      if (buffer_includes_halos) then
        i0 = isc - isd + 1
      endif
      i1 = i0 + xc_size - 1
    case (2)
      ! X is the second dimension
      if (buffer_includes_halos) then
        j0 = isc - isd + 1
      endif
      j1 = j0 + xc_size - 1
    case (3)
      ! X is the third dimension
      if (buffer_includes_halos) then
        k0 = isc - isd + 1
      endif
      k1 = k0 + xc_size - 1
  end select

  select case (ydim_index)
    case (1)
      ! y is the first dimension
      if (buffer_includes_halos) then
        i0 = jsc - jsd + 1
      endif
      i1 = i0 + yc_size - 1
    case (2)
      ! y is the second dimension
      if (buffer_includes_halos) then
        j0 = jsc - jsd + 1
      endif
      j1 = j0 + yc_size - 1
    case (3)
      ! y is the third dimension
      if (buffer_includes_halos) then
        k0 = jsc - jsd + 1
      endif
      k1 = k0 + yc_size - 1
  end select

  varid = get_variable_id(fileobj%ncid, trim(variable_name), &
    msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))

  select type(vdata)
    type is (integer(kind=i4_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1, k0:k1), &
        start=c, count=e)
    type is (integer(kind=i8_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1, k0:k1), &
        start=c, count=e)
    type is (real(kind=r4_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1, k0:k1), &
        start=c, count=e)
    type is (real(kind=r8_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1, k0:k1), &
        start=c, count=e)
  end select

  call check_netcdf_code(err, "Failed to write variable: "//trim(variable_name))
end subroutine netcdf_mpi_write_3d

!> @brief Write 4D data using NetCDF MPI
subroutine netcdf_mpi_write_4d(fileobj, variable_name, vdata, unlim_dim_level, corner, edge_lengths)
  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(in) :: vdata !< Data that will be written out to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited dimension.
  integer, dimension(4), intent(in), optional :: corner !< Array of starting indices describing where the data
                                                        !! will be written to.
  integer, dimension(4), intent(in), optional :: edge_lengths !< The number of elements that will be written
                                                              !! in each dimension.

  integer :: xdim_index
  integer :: ydim_index
  integer :: xpos
  integer :: ypos

  integer :: isd !< Starting "x" index of the data domain
  integer :: isc !< Starting "x" index of the compute domain
  integer :: jsd !< Starting "y" index of the data domain
  integer :: jsc !< Starting "y" index of the compute domain
  integer :: xc_size !< Size of the compute domain
  integer :: yc_size !< Size of the compute domain
  logical :: buffer_includes_halos !< True if "vdata" is the size of the data domain
                                   !! (i.e it includes halos (which are not written))

  integer, dimension(5) :: c !< Indices of the corners for each dimension
  integer, dimension(5) :: e !< Size of the data for each dimension
  integer :: varid
  integer :: unlim_dim_index

  integer :: i0, i1
  integer :: j0, j1
  integer :: k0, k1
  integer :: l0, l1
  integer :: err

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., xdim_index, ydim_index, xpos, ypos)) then
    call compressed_write(fileobj, variable_name, vdata, unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)
    return
  endif

  !! Get some more info about the variable
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, &
                      buffer_includes_halos, msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))

  c = 1
  e = [shape(vdata), 1]

  c(xdim_index) = isc
  c(ydim_index) = jsc
  e(xdim_index) = xc_size
  e(ydim_index) = yc_size

  if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, broadcast=.false.)
      if (unlim_dim_index .ne. 5) then
        call error("unlimited dimension must be the slowest varying dimension in variable: "//trim(variable_name))
      endif
      c(unlim_dim_index) = unlim_dim_level
  endif

  i0 = 1
  j0 = 1
  k0 = 1
  l0 = 1

  i1 = size(vdata, 1)
  j1 = size(vdata, 2)
  k1 = size(vdata, 3)
  l1 = size(vdata, 4)

  select case (xdim_index)
    case (1)
      ! X is the first dimension
      if (buffer_includes_halos) then
        i0 = isc - isd + 1
      endif
      i1 = i0 + xc_size - 1
    case (2)
      ! X is the second dimension
      if (buffer_includes_halos) then
        j0 = isc - isd + 1
      endif
      j1 = j0 + xc_size - 1
    case (3)
      ! X is the third dimension
      if (buffer_includes_halos) then
        k0 = isc - isd + 1
      endif
      k1 = k0 + xc_size - 1
    case (4)
      ! X is the fourth dimension
      if (buffer_includes_halos) then
        l0 = isc - isd + 1
      endif
      l1 = l0 + xc_size - 1
  end select

  select case (ydim_index)
    case (1)
      ! y is the first dimension
      if (buffer_includes_halos) then
        i0 = jsc - jsd + 1
      endif
      i1 = i0 + yc_size - 1
    case (2)
      ! y is the second dimension
      if (buffer_includes_halos) then
        j0 = jsc - jsd + 1
      endif
      j1 = j0 + yc_size - 1
    case (3)
      ! y is the third dimension
      if (buffer_includes_halos) then
        k0 = jsc - jsd + 1
      endif
      k1 = k0 + yc_size - 1
    case (4)
      ! y is the fourth dimension
      if (buffer_includes_halos) then
        l0 = jsc - jsd + 1
      endif
      l1 = l0 + yc_size - 1
  end select

  varid = get_variable_id(fileobj%ncid, trim(variable_name), &
    msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))

  select type(vdata)
    type is (integer(kind=i4_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1, k0:k1, l0:l1), &
        start=c, count=e)
    type is (integer(kind=i8_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1, k0:k1, l0:l1), &
        start=c, count=e)
    type is (real(kind=r4_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1, k0:k1, l0:l1), &
        start=c, count=e)
    type is (real(kind=r8_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1, k0:k1, l0:l1), &
        start=c, count=e)
  end select

  call check_netcdf_code(err, "Failed to write variable: "//trim(variable_name))
end subroutine netcdf_mpi_write_4d

!> @brief Write 5D data using NetCDF MPI
subroutine netcdf_mpi_write_5d(fileobj, variable_name, vdata, unlim_dim_level, corner, edge_lengths)
  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(in) :: vdata !< Data that will be written out to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited dimension.
  integer, dimension(5), intent(in), optional :: corner !< Array of starting indices describing where the data
                                                        !! will be written to.
  integer, dimension(5), intent(in), optional :: edge_lengths !< The number of elements that will be written
                                                              !! in each dimension.

  integer :: xdim_index
  integer :: ydim_index
  integer :: xpos
  integer :: ypos

  integer :: isd !< Starting "x" index of the data domain
  integer :: isc !< Starting "x" index of the compute domain
  integer :: jsd !< Starting "y" index of the data domain
  integer :: jsc !< Starting "y" index of the compute domain
  integer :: xc_size !< Size of the compute domain
  integer :: yc_size !< Size of the compute domain
  logical :: buffer_includes_halos !< True if "vdata" is the size of the data domain
                                   !! (i.e it includes halos (which are not written))

  integer, dimension(6) :: c !< Indices of the corners for each dimension
  integer, dimension(6) :: e !< Size of the data for each dimension
  integer :: varid
  integer :: unlim_dim_index

  integer :: i0, i1
  integer :: j0, j1
  integer :: k0, k1
  integer :: l0, l1
  integer :: m0, m1
  integer :: err

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., xdim_index, ydim_index, xpos, ypos)) then
    call compressed_write(fileobj, variable_name, vdata, unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths)
    return
  endif

  !! Get some more info about the variable
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, &
                      buffer_includes_halos, msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))

  c = 1
  e = [shape(vdata), 1]

  c(xdim_index) = isc
  c(ydim_index) = jsc
  e(xdim_index) = xc_size
  e(ydim_index) = yc_size

  if (present(unlim_dim_level)) then
      unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, variable_name, broadcast=.false.)
      if (unlim_dim_index .ne. 6) then
        call error("unlimited dimension must be the slowest varying dimension in variable: "//trim(variable_name))
      endif
      c(unlim_dim_index) = unlim_dim_level
  endif

  i0 = 1
  j0 = 1
  k0 = 1
  l0 = 1
  m0 = 1

  i1 = size(vdata, 1)
  j1 = size(vdata, 2)
  k1 = size(vdata, 3)
  l1 = size(vdata, 4)
  m1 = size(vdata, 5)

  select case (xdim_index)
    case (1)
      ! X is the first dimension
      if (buffer_includes_halos) then
        i0 = isc - isd + 1
      endif
      i1 = i0 + xc_size - 1
    case (2)
      ! X is the second dimension
      if (buffer_includes_halos) then
        j0 = isc - isd + 1
      endif
      j1 = j0 + xc_size - 1
    case (3)
      ! X is the third dimension
      if (buffer_includes_halos) then
        k0 = isc - isd + 1
      endif
      k1 = k0 + xc_size - 1
    case (4)
      ! X is the fourth dimension
      if (buffer_includes_halos) then
        l0 = isc - isd + 1
      endif
      l1 = l0 + xc_size - 1
    case (5)
      ! X is the fifth dimension
      if (buffer_includes_halos) then
        m0 = isc - isd + 1
      endif
      m1 = m0 + xc_size - 1
  end select

  select case (ydim_index)
    case (1)
      ! y is the first dimension
      if (buffer_includes_halos) then
        i0 = jsc - jsd + 1
      endif
      i1 = i0 + yc_size - 1
    case (2)
      ! y is the second dimension
      if (buffer_includes_halos) then
        j0 = jsc - jsd + 1
      endif
      j1 = j0 + yc_size - 1
    case (3)
      ! y is the third dimension
      if (buffer_includes_halos) then
        k0 = jsc - jsd + 1
      endif
      k1 = k0 + yc_size - 1
    case (4)
      ! y is the fourth dimension
      if (buffer_includes_halos) then
        l0 = jsc - jsd + 1
      endif
      l1 = l0 + yc_size - 1
    case (5)
      ! y is the fifth dimension
      if (buffer_includes_halos) then
        m0 = jsc - jsd + 1
      endif
      m1 = m0 + yc_size - 1
  end select

  varid = get_variable_id(fileobj%ncid, trim(variable_name), &
    msg="file:"//trim(fileobj%path)//" and variable:"//trim(variable_name))

  select type(vdata)
    type is (integer(kind=i4_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1, k0:k1, l0:l1, m0:m1), &
        start=c, count=e)
    type is (integer(kind=i8_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1, k0:k1, l0:l1, m0:m1), &
        start=c, count=e)
    type is (real(kind=r4_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1, k0:k1, l0:l1, m0:m1), &
        start=c, count=e)
    type is (real(kind=r8_kind))
      err = nf90_put_var(fileobj%ncid, varid, &
        vdata(i0:i1, j0:j1, k0:k1, l0:l1, m0:m1), &
        start=c, count=e)
  end select

  call check_netcdf_code(err, "Failed to write variable: "//trim(variable_name))
end subroutine netcdf_mpi_write_5d
!> @}
