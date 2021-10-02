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
!> @file
!> @ingroup fms_netcdf_domain_io_mod

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
!!        decomposed".
subroutine domain_write_2d(fileobj, variable_name, vdata, unlim_dim_level, &
                           corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(in) :: vdata !< Data that will
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

  integer(kind=i4_kind), dimension(:,:), allocatable :: buf_i4_kind
  integer(kind=i8_kind), dimension(:,:), allocatable :: buf_i8_kind
  real(kind=r4_kind), dimension(:,:), allocatable :: buf_r4_kind
  real(kind=r8_kind), dimension(:,:), allocatable :: buf_r8_kind
  logical :: buffer_includes_halos
  integer, dimension(2) :: c
  integer, dimension(2) :: e
  integer(kind=i4_kind), dimension(:,:), allocatable :: global_buf_i4_kind
  integer(kind=i8_kind), dimension(:,:), allocatable :: global_buf_i8_kind
  real(kind=r4_kind), dimension(:,:), allocatable :: global_buf_r4_kind
  real(kind=r8_kind), dimension(:,:), allocatable :: global_buf_r8_kind
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
  real(kind=i4_kind) :: fill_i4_kind !< Fill value of a i4_kind variable
  real(kind=i8_kind) :: fill_i8_kind !< Fill value of a i8_kind variable
  real(kind=r4_kind) :: fill_r4_kind !< Fill value of a r4_kind variable
  real(kind=r8_kind) :: fill_r8_kind !< Fill value of a r8_kind variable
  integer :: xgmax !< Ending x index of the global io domain
  integer :: xgmin !< Starting x index of the global io domain
  integer :: ygmax !< Ending y index of the global io domain
  integer :: ygmin !< Ending y index of the global io domain

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
      call error("unsupported variable type: domain_write_2d: file: "//trim(fileobj%path)//" variable:"//trim(variable_name))
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
        call error("unsupported variable type: domain_write_2d: file: "//trim(fileobj%path)//" variable:"//trim(variable_name))
    end select
  endif
end subroutine domain_write_2d


!> @brief Gather "compute" domain data on the I/O root rank and then have
!!        the I/O root write out the data that spans the "global" domain.
!!        This routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_write_3d(fileobj, variable_name, vdata, unlim_dim_level, &
                           corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(in) :: vdata !< Data that will
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

  integer(kind=i4_kind), dimension(:,:,:), allocatable :: buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:), allocatable :: buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:), allocatable :: buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:), allocatable :: buf_r8_kind
  logical :: buffer_includes_halos
  integer, dimension(3) :: c
  integer, dimension(3) :: e
  integer(kind=i4_kind), dimension(:,:,:), allocatable :: global_buf_i4_kind
  integer(kind=i8_kind), dimension(:,:,:), allocatable :: global_buf_i8_kind
  real(kind=r4_kind), dimension(:,:,:), allocatable :: global_buf_r4_kind
  real(kind=r8_kind), dimension(:,:,:), allocatable :: global_buf_r8_kind
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
  real(kind=i4_kind) :: fill_i4_kind !< Fill value of a i4_kind variable
  real(kind=i8_kind) :: fill_i8_kind !< Fill value of a i8_kind variable
  real(kind=r4_kind) :: fill_r4_kind !< Fill value of a r4_kind variable
  real(kind=r8_kind) :: fill_r8_kind !< Fill value of a r8_kind variable
  integer :: xgmax !< Ending x index of the global io domain
  integer :: xgmin !< Starting x index of the global io domain
  integer :: ygmax !< Ending y index of the global io domain
  integer :: ygmin !< Ending y index of the global io domain

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
        call error("unsupported variable type: domain_write_3d: file: "//trim(fileobj%path)//" variable:"//trim(variable_name))
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
        call error("unsupported variable type: domain_write_3d: file: "//trim(fileobj%path)//" variable:"//trim(variable_name))
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
  real(kind=i4_kind) :: fill_i4_kind !< Fill value of a i4_kind variable
  real(kind=i8_kind) :: fill_i8_kind !< Fill value of a i8_kind variable
  real(kind=r4_kind) :: fill_r4_kind !< Fill value of a r4_kind variable
  real(kind=r8_kind) :: fill_r8_kind !< Fill value of a r8_kind variable
  integer :: xgmax !< Ending x index of the global io domain
  integer :: xgmin !< Starting x index of the global io domain
  integer :: ygmax !< Ending y index of the global io domain
  integer :: ygmin !< Ending y index of the global io domain

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
        call error("unsupported variable type: domain_write_4d: file: "//trim(fileobj%path)//" variable:"//trim(variable_name))
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
        call error("unsupported variable type: domain_write_4d: file: "//trim(fileobj%path)//" variable:"//trim(variable_name))
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
  real(kind=i4_kind) :: fill_i4_kind !< Fill value of a i4_kind variable
  real(kind=i8_kind) :: fill_i8_kind !< Fill value of a i8_kind variable
  real(kind=r4_kind) :: fill_r4_kind !< Fill value of a r4_kind variable
  real(kind=r8_kind) :: fill_r8_kind !< Fill value of a r8_kind variable
  integer :: xgmax !< Ending x index of the global io domain
  integer :: xgmin !< Starting x index of the global io domain
  integer :: ygmax !< Ending y index of the global io domain
  integer :: ygmin !< Ending y index of the global io domain

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
        call error("unsupported variable type: domain_write_5d: file: "//trim(fileobj%path)//" variable:"//trim(variable_name))
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
        call error("unsupported variable type: domain_write_5d: file: "//trim(fileobj%path)//" variable:"//trim(variable_name))
    end select
  endif
end subroutine domain_write_5d
