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

!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_read_0d(fileobj, variable_name, vdata, unlim_dim_level, corner)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), intent(inout) :: vdata !< Data that will
                                   !! be written out
                                   !! to the netcdf file.
  integer, intent(in), optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.
  integer, intent(in), optional :: corner !< Array of starting
                                          !! indices describing
                                          !! where the data
                                          !! will be written to.

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true.)) then
    call netcdf_read_data(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          broadcast=.true.)
    return
  else
    call error("this branch should never be reached.")
  endif
end subroutine domain_read_0d


!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_read_1d(fileobj, variable_name, vdata, unlim_dim_level, &
                          corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:), intent(inout) :: vdata !< Data that will
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

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true.)) then
    call netcdf_read_data(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths, broadcast=.true.)
    return
  else
    call error("this branch should never be reached.")
  endif
end subroutine domain_read_1d


!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_read_2d(fileobj, variable_name, vdata, unlim_dim_level, &
                          corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:), intent(inout) :: vdata !< Data that will
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

  integer :: xdim_index
  integer :: ydim_index
  type(domain2d), pointer :: io_domain
  integer :: xpos
  integer :: ypos
  integer :: i
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_jsc
  integer, dimension(:), allocatable :: pe_jcsize
  integer, dimension(2) :: c
  integer, dimension(2) :: e
  integer(kind=int32), dimension(:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:), allocatable :: buf_real64
  logical :: buffer_includes_halos
  integer :: xgmin !< Starting x index of global io domain
  integer :: ygmin !< Starting y index of global io domain

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index, xpos, ypos)) then
    call netcdf_read_data(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths, broadcast=.true.)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, buffer_includes_halos)
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root reads in the data and scatters it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xsize=pe_icsize, position=xpos)
    call mpp_get_compute_domains(io_domain, ybegin=pe_jsc, ysize=pe_jcsize, position=ypos)
    call mpp_get_global_domain(io_domain, xbegin=xgmin, position=xpos)
    call mpp_get_global_domain(io_domain, ybegin=ygmin, position=ypos)
    do i = 1, size(fileobj%pelist)
      c(xdim_index) = pe_isc(i)
      c(ydim_index) = pe_jsc(i)
      if (fileobj%adjust_indices) then
        c(xdim_index) = c(xdim_index) - xgmin + 1
        c(ydim_index) = c(ydim_index) - ygmin + 1
      endif
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=int32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int32, e)
          call netcdf_read_data(fileobj, variable_name, buf_int32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_int32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int32, size(buf_int32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int32)
        type is (integer(kind=int64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int64, e)
          call netcdf_read_data(fileobj, variable_name, buf_int64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_int64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int64, size(buf_int64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int64)
        type is (real(kind=real32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real32, e)
          call netcdf_read_data(fileobj, variable_name, buf_real32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_real32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real32, size(buf_real32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real32)
        type is (real(kind=real64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real64, e)
          call netcdf_read_data(fileobj, variable_name, buf_real64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_real64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real64, size(buf_real64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real64)
        class default
          call error("domain_read_2d: Unsupported type for variable: "//variable_name//" in file: "//trim(fileobj%path)//"")
      end select
    enddo
    deallocate(pe_isc)
    deallocate(pe_icsize)
    deallocate(pe_jsc)
    deallocate(pe_jcsize)
  else
    if (buffer_includes_halos) then
      c(xdim_index) = isc - isd + 1
      c(ydim_index) = jsc - jsd + 1
    endif
    e(xdim_index) = xc_size
    e(ydim_index) = yc_size
    select type(vdata)
      type is (integer(kind=int32))
        call allocate_array(buf_int32, e)
        call mpp_recv(buf_int32, size(buf_int32), fileobj%io_root, block=.true.)
        call put_array_section(buf_int32, vdata, c, e)
        deallocate(buf_int32)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, e)
        call mpp_recv(buf_int64, size(buf_int64), fileobj%io_root, block=.true.)
        call put_array_section(buf_int64, vdata, c, e)
        deallocate(buf_int64)
      type is (real(kind=real32))
        call allocate_array(buf_real32, e)
        call mpp_recv(buf_real32, size(buf_real32), fileobj%io_root, block=.true.)
        call put_array_section(buf_real32, vdata, c, e)
        deallocate(buf_real32)
      type is (real(kind=real64))
        call allocate_array(buf_real64, e)
        call mpp_recv(buf_real64, size(buf_real64), fileobj%io_root, block=.true.)
        call put_array_section(buf_real64, vdata, c, e)
        deallocate(buf_real64)
      class default
          call error("domain_read_2d: Unsupported type for variable: "//variable_name//" in file: "//trim(fileobj%path)//"")
    end select
  endif
end subroutine domain_read_2d


!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_read_3d(fileobj, variable_name, vdata, unlim_dim_level, &
                          corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:), intent(inout) :: vdata !< Data that will
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

  integer :: xdim_index
  integer :: ydim_index
  type(domain2d), pointer :: io_domain
  integer :: xpos
  integer :: ypos
  integer :: i
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_jsc
  integer, dimension(:), allocatable :: pe_jcsize
  integer, dimension(3) :: c
  integer, dimension(3) :: e
  integer(kind=int32), dimension(:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:), allocatable :: buf_real64
  logical :: buffer_includes_halos
  integer :: xgmin !< Starting x index of global io domain
  integer :: ygmin !< Starting y index of global io domain

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index, xpos, ypos)) then
    call netcdf_read_data(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths, broadcast=.true.)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, buffer_includes_halos)
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root reads in the data and scatters it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xsize=pe_icsize, position=xpos)
    call mpp_get_compute_domains(io_domain, ybegin=pe_jsc, ysize=pe_jcsize, position=ypos)
    call mpp_get_global_domain(io_domain, xbegin=xgmin, position=xpos)
    call mpp_get_global_domain(io_domain, ybegin=ygmin, position=ypos)
    do i = 1, size(fileobj%pelist)
      c(xdim_index) = pe_isc(i)
      c(ydim_index) = pe_jsc(i)
      if (fileobj%adjust_indices) then
        c(xdim_index) = c(xdim_index) - xgmin + 1
        c(ydim_index) = c(ydim_index) - ygmin + 1
      endif
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=int32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int32, e)
          call netcdf_read_data(fileobj, variable_name, buf_int32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_int32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int32, size(buf_int32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int32)
        type is (integer(kind=int64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int64, e)
          call netcdf_read_data(fileobj, variable_name, buf_int64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_int64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int64, size(buf_int64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int64)
        type is (real(kind=real32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real32, e)
          call netcdf_read_data(fileobj, variable_name, buf_real32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_real32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real32, size(buf_real32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real32)
        type is (real(kind=real64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real64, e)
          call netcdf_read_data(fileobj, variable_name, buf_real64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_real64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real64, size(buf_real64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real64)
        class default
          call error("domain_read_3d: Unsupported type for variable: "//variable_name//" in file: "//trim(fileobj%path)//"")
      end select
    enddo
    deallocate(pe_isc)
    deallocate(pe_icsize)
    deallocate(pe_jsc)
    deallocate(pe_jcsize)
  else
    if (buffer_includes_halos) then
      c(xdim_index) = isc - isd + 1
      c(ydim_index) = jsc - jsd + 1
    endif
    e(xdim_index) = xc_size
    e(ydim_index) = yc_size
    select type(vdata)
      type is (integer(kind=int32))
        call allocate_array(buf_int32, e)
        call mpp_recv(buf_int32, size(buf_int32), fileobj%io_root, block=.true.)
        call put_array_section(buf_int32, vdata, c, e)
        deallocate(buf_int32)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, e)
        call mpp_recv(buf_int64, size(buf_int64), fileobj%io_root, block=.true.)
        call put_array_section(buf_int64, vdata, c, e)
        deallocate(buf_int64)
      type is (real(kind=real32))
        call allocate_array(buf_real32, e)
        call mpp_recv(buf_real32, size(buf_real32), fileobj%io_root, block=.true.)
        call put_array_section(buf_real32, vdata, c, e)
        deallocate(buf_real32)
      type is (real(kind=real64))
        call allocate_array(buf_real64, e)
        call mpp_recv(buf_real64, size(buf_real64), fileobj%io_root, block=.true.)
        call put_array_section(buf_real64, vdata, c, e)
        deallocate(buf_real64)
      class default
          call error("domain_read_3d: Unsupported type for variable: "//variable_name//" in file: "//trim(fileobj%path)//"")
    end select
  endif
end subroutine domain_read_3d


!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_read_4d(fileobj, variable_name, vdata, unlim_dim_level, &
                          corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:), intent(inout) :: vdata !< Data that will
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

  integer :: xdim_index
  integer :: ydim_index
  type(domain2d), pointer :: io_domain
  integer :: xpos
  integer :: ypos
  integer :: i
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_jsc
  integer, dimension(:), allocatable :: pe_jcsize
  integer, dimension(4) :: c
  integer, dimension(4) :: e
  integer(kind=int32), dimension(:,:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:,:), allocatable :: buf_real64
  logical :: buffer_includes_halos
  integer :: xgmin !< Starting x index of global io domain
  integer :: ygmin !< Starting y index of global io domain

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index, xpos, ypos)) then
    call netcdf_read_data(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths, broadcast=.true.)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, buffer_includes_halos)
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root reads in the data and scatters it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xsize=pe_icsize, position=xpos)
    call mpp_get_compute_domains(io_domain, ybegin=pe_jsc, ysize=pe_jcsize, position=ypos)
    call mpp_get_global_domain(io_domain, xbegin=xgmin, position=xpos)
    call mpp_get_global_domain(io_domain, ybegin=ygmin, position=ypos)
    do i = 1, size(fileobj%pelist)
      c(xdim_index) = pe_isc(i)
      c(ydim_index) = pe_jsc(i)
      if (fileobj%adjust_indices) then
        c(xdim_index) = c(xdim_index) - xgmin + 1
        c(ydim_index) = c(ydim_index) - ygmin + 1
      endif
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=int32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int32, e)
          call netcdf_read_data(fileobj, variable_name, buf_int32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_int32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int32, size(buf_int32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int32)
        type is (integer(kind=int64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int64, e)
          call netcdf_read_data(fileobj, variable_name, buf_int64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_int64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int64, size(buf_int64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int64)
        type is (real(kind=real32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real32, e)
          call netcdf_read_data(fileobj, variable_name, buf_real32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_real32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real32, size(buf_real32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real32)
        type is (real(kind=real64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real64, e)
          call netcdf_read_data(fileobj, variable_name, buf_real64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_real64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real64, size(buf_real64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real64)
        class default
          call error("domain_read_4d: Unsupported type for variable: "//variable_name//" in file: "//trim(fileobj%path)//"")
      end select
    enddo
    deallocate(pe_isc)
    deallocate(pe_icsize)
    deallocate(pe_jsc)
    deallocate(pe_jcsize)
  else
    if (buffer_includes_halos) then
      c(xdim_index) = isc - isd + 1
      c(ydim_index) = jsc - jsd + 1
    endif
    e(xdim_index) = xc_size
    e(ydim_index) = yc_size
    select type(vdata)
      type is (integer(kind=int32))
        call allocate_array(buf_int32, e)
        call mpp_recv(buf_int32, size(buf_int32), fileobj%io_root, block=.true.)
        call put_array_section(buf_int32, vdata, c, e)
        deallocate(buf_int32)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, e)
        call mpp_recv(buf_int64, size(buf_int64), fileobj%io_root, block=.true.)
        call put_array_section(buf_int64, vdata, c, e)
        deallocate(buf_int64)
      type is (real(kind=real32))
        call allocate_array(buf_real32, e)
        call mpp_recv(buf_real32, size(buf_real32), fileobj%io_root, block=.true.)
        call put_array_section(buf_real32, vdata, c, e)
        deallocate(buf_real32)
      type is (real(kind=real64))
        call allocate_array(buf_real64, e)
        call mpp_recv(buf_real64, size(buf_real64), fileobj%io_root, block=.true.)
        call put_array_section(buf_real64, vdata, c, e)
        deallocate(buf_real64)
      class default
          call error("domain_read_4d: Unsupported type for variable: "//variable_name//" in file: "//trim(fileobj%path)//"")
    end select
  endif
end subroutine domain_read_4d


!> @brief I/O domain root reads in  a domain decomposed variable at a
!!        specific unlimited dimension level and scatters the data to the
!!        rest of the ranks using its I/O compute domain indices. This
!!        routine may only be used with variables that are "domain
!!        decomposed".
subroutine domain_read_5d(fileobj, variable_name, vdata, unlim_dim_level, &
                          corner, edge_lengths)

  type(FmsNetcdfDomainFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  class(*), dimension(:,:,:,:,:), intent(inout) :: vdata !< Data that will
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

  integer :: xdim_index
  integer :: ydim_index
  type(domain2d), pointer :: io_domain
  integer :: xpos
  integer :: ypos
  integer :: i
  integer :: isd
  integer :: isc
  integer :: xc_size
  integer :: jsd
  integer :: jsc
  integer :: yc_size
  integer, dimension(:), allocatable :: pe_isc
  integer, dimension(:), allocatable :: pe_icsize
  integer, dimension(:), allocatable :: pe_jsc
  integer, dimension(:), allocatable :: pe_jcsize
  integer, dimension(5) :: c
  integer, dimension(5) :: e
  integer(kind=int32), dimension(:,:,:,:,:), allocatable :: buf_int32
  integer(kind=int64), dimension(:,:,:,:,:), allocatable :: buf_int64
  real(kind=real32), dimension(:,:,:,:,:), allocatable :: buf_real32
  real(kind=real64), dimension(:,:,:,:,:), allocatable :: buf_real64
  logical :: buffer_includes_halos
  integer :: xgmin !< Starting x index of global io domain
  integer :: ygmin !< Starting y index of global io domain

  if (.not. is_variable_domain_decomposed(fileobj, variable_name, .true., &
                                          xdim_index, ydim_index, xpos, ypos)) then
    call netcdf_read_data(fileobj, variable_name, vdata, &
                          unlim_dim_level=unlim_dim_level, corner=corner, &
                          edge_lengths=edge_lengths, broadcast=.true.)
    return
  endif
  io_domain => mpp_get_io_domain(fileobj%domain)
  call domain_offsets(size(vdata, xdim_index), size(vdata, ydim_index), fileobj%domain, &
                      xpos, ypos, isd, isc, xc_size, jsd, jsc, yc_size, buffer_includes_halos)
  c(:) = 1
  e(:) = shape(vdata)

  !I/O root reads in the data and scatters it.
  if (fileobj%is_root) then
    allocate(pe_isc(size(fileobj%pelist)))
    allocate(pe_icsize(size(fileobj%pelist)))
    allocate(pe_jsc(size(fileobj%pelist)))
    allocate(pe_jcsize(size(fileobj%pelist)))
    call mpp_get_compute_domains(io_domain, xbegin=pe_isc, xsize=pe_icsize, position=xpos)
    call mpp_get_compute_domains(io_domain, ybegin=pe_jsc, ysize=pe_jcsize, position=ypos)
    call mpp_get_global_domain(io_domain, xbegin=xgmin, position=xpos)
    call mpp_get_global_domain(io_domain, ybegin=ygmin, position=ypos)
    do i = 1, size(fileobj%pelist)
      !Calculate the indices of the domain-decomposed chunk relative to its position in the file.
      c(xdim_index) = pe_isc(i)
      c(ydim_index) = pe_jsc(i)
      if (fileobj%adjust_indices) then
        c(xdim_index) = c(xdim_index) - xgmin + 1
        c(ydim_index) = c(ydim_index) - ygmin + 1
      endif
      e(xdim_index) = pe_icsize(i)
      e(ydim_index) = pe_jcsize(i)
      select type(vdata)
        type is (integer(kind=int32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int32, e)
          call netcdf_read_data(fileobj, variable_name, buf_int32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.  Re-adjust the indicies relative
            !to the input buffer vdata.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_int32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int32, size(buf_int32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int32)
        type is (integer(kind=int64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_int64, e)
          call netcdf_read_data(fileobj, variable_name, buf_int64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_int64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_int64, size(buf_int64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_int64)
        type is (real(kind=real32))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real32, e)
          call netcdf_read_data(fileobj, variable_name, buf_real32, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_real32, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real32, size(buf_real32), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real32)
        type is (real(kind=real64))
          !Read in the data for fileobj%pelist(i)'s portion of the compute domain.
          call allocate_array(buf_real64, e)
          call netcdf_read_data(fileobj, variable_name, buf_real64, &
                                unlim_dim_level=unlim_dim_level, &
                                corner=c, edge_lengths=e, broadcast=.false.)
          if (i .eq. 1) then
            !Root rank stores data directly.
            if (buffer_includes_halos) then
              !Adjust if the input buffer has room for halos.
              c(xdim_index) = isc - isd + 1
              c(ydim_index) = jsc - jsd + 1
            else
              c(xdim_index) = 1
              c(ydim_index) = 1
            endif
            call put_array_section(buf_real64, vdata, c, e)
          else
            !Send data to non-root ranks.
            call mpp_send(buf_real64, size(buf_real64), fileobj%pelist(i))
            call mpp_sync_self(check=EVENT_SEND)
          endif
          deallocate(buf_real64)
        class default
          call error("domain_read_5d: Unsupported type for variable: "//variable_name//" in file: "//trim(fileobj%path)//"")
      end select
    enddo
    deallocate(pe_isc)
    deallocate(pe_icsize)
    deallocate(pe_jsc)
    deallocate(pe_jcsize)
  else
    if (buffer_includes_halos) then
      c(xdim_index) = isc - isd + 1
      c(ydim_index) = jsc - jsd + 1
    endif
    e(xdim_index) = xc_size
    e(ydim_index) = yc_size
    select type(vdata)
      type is (integer(kind=int32))
        call allocate_array(buf_int32, e)
        call mpp_recv(buf_int32, size(buf_int32), fileobj%io_root, block=.true.)
        call put_array_section(buf_int32, vdata, c, e)
        deallocate(buf_int32)
      type is (integer(kind=int64))
        call allocate_array(buf_int64, e)
        call mpp_recv(buf_int64, size(buf_int64), fileobj%io_root, block=.true.)
        call put_array_section(buf_int64, vdata, c, e)
        deallocate(buf_int64)
      type is (real(kind=real32))
        call allocate_array(buf_real32, e)
        call mpp_recv(buf_real32, size(buf_real32), fileobj%io_root, block=.true.)
        call put_array_section(buf_real32, vdata, c, e)
        deallocate(buf_real32)
      type is (real(kind=real64))
        call allocate_array(buf_real64, e)
        call mpp_recv(buf_real64, size(buf_real64), fileobj%io_root, block=.true.)
        call put_array_section(buf_real64, vdata, c, e)
        deallocate(buf_real64)
      class default
          call error("domain_read_5d: Unsupported type for variable: "//variable_name//" in file: "//trim(fileobj%path)//"")
    end select
  endif
end subroutine domain_read_5d
