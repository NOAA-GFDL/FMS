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
!> @ingroup fms2_io

!> @brief gathers the 2d vdata from all of the relevant pes into the root_pe and saves it to a
!! buffer.
subroutine gather_data_bc_2d(fileobj, vdata, bc_info)
  class(FmsNetcdfFile_t), intent(inout) :: fileobj !< Fms2io netcdf fileobj
  class(*), dimension(:,:), intent(in) :: vdata !< Data to be gather
  type(bc_information), intent(inout) :: bc_info !< information about the boundary condition variable

  integer :: i_glob !< Size of the global domain in x
  integer :: j_glob !< Size of the global domain in y

  integer :: isc, iec, jsc, jec !< current PE's indices
  integer :: i1, i2, j1, j2 !< current PE's indices relative to the global domain
  integer :: i_add, j_add !< current PE's shift

  real(kind=r4_kind), dimension(:,:), allocatable, target :: global_buf_r4_kind !< buffer with a data gathered
  real(kind=r4_kind), dimension(:,:), allocatable, target :: local_buf_r4_kind !< current PEs data
  real(kind=r8_kind), dimension(:,:), allocatable, target :: global_buf_r8_kind !< buffer with a data gathered
  real(kind=r8_kind), dimension(:,:), allocatable, target :: local_buf_r8_kind !< current PEs data

  integer(kind=i8_kind) :: chksum_val !< Checksum value calculated from the gathered data
  character(len=32) :: chksum !< Cheksum value converted to a string

  !> Set the indices
  isc = bc_info%indices(1)
  iec = bc_info%indices(2)
  jsc = bc_info%indices(3)
  jec = bc_info%indices(4)

  !> This is the section of the PE that will actually be added to the global_buffer
  i1 = 1 + bc_info%x_halo
  i2 = i1 + (iec-isc)
  j1 = 1 + bc_info%y_halo
  j2 = j1 + (jec-jsc)

  !> Set up index shifts for global array
  i_add = bc_info%ishift
  j_add = bc_info%jshift

  !> Only root allocates a global_buffer that will be written
  if (fileobj%is_root) then
        i_glob = bc_info%global_size(1)
        j_glob = bc_info%global_size(2)
  endif

  !> Gather the data and calculate the checksum for the resulting array.
  select type(vdata)
    type is (real(kind=r4_kind))
       !> If the fileobj's root does not have data for this variable
       if (fileobj%is_root .and. .not. bc_info%data_on_file_root) then
           !> Allocate global_buf_r4_kind to be one size bigger, global_buf_r4_kind(i_glob+1,,:) is just dummy data
          allocate(global_buf_r4_kind(i_glob+1, j_glob))
          !> Allocate a temp local buffer to the fileobj's root. This is needed because the data needs to be send
          !! to the fileobj's root, but because fileobj's root doesn't have any data, we just create dummy data and
          !! not write it later
          allocate(local_buf_r4_kind(1,1))
          local_buf_r4_kind = 0.
          isc = 1+i_glob; i_add=0; iec=1+i_glob; jsc=j_glob; j_add=0; jec=j_glob
          i1=1; i2=1; j1=1; j2=1
       else
          !! In this case there is data in fileobj's root, so there is no need for the dummy data
          if(fileobj%is_root) allocate(global_buf_r4_kind(i_glob, j_glob))
          allocate(local_buf_r4_kind(size(vdata,1), size(vdata,2)))
          local_buf_r4_kind = vdata
       endif

       call mpp_gather(isc+i_add, iec+i_add, jsc+j_add, jec+j_add, bc_info%pelist, &
                          local_buf_r4_kind(i1:i2,j1:j2), &
                          global_buf_r4_kind, fileobj%is_root)

        deallocate(local_buf_r4_kind)
       !> If you are on fileobj's root, calculate the checksum and save the gathered data in a buffer
       if (fileobj%is_root) then
           chksum_val = mpp_chksum(global_buf_r4_kind(1:i_glob,1:j_glob), (/mpp_pe()/))
           allocate(bc_info%globaldata2d_r4(i_glob, j_glob))
           bc_info%globaldata2d_r4=global_buf_r4_kind(1:i_glob,1:j_glob)
           deallocate(global_buf_r4_kind)
       endif
    type is (real(kind=r8_kind))
       !> If the fileobj's root does not have data for this variable
       if (fileobj%is_root .and. .not. bc_info%data_on_file_root) then
          !> Allocate global_buf_r8_kind to be one size bigger, global_buf_r8_kind(i_glob+1,,:) is just dummy data
          allocate(global_buf_r8_kind(i_glob+1, j_glob))
          !> Allocate a temp local buffer to the fileobj's root. This is needed because the data needs to be send
          !! to the fileobj's root, but because fileobj's root doesn't have any data, we just create dummy data and
          !! not write it later
          allocate(local_buf_r8_kind(1,1))
          local_buf_r8_kind = 0.
          isc = 1+i_glob; i_add=0; iec=1+i_glob; jsc=j_glob; j_add=0; jec=j_glob
          i1=1; i2=1; j1=1; j2=1
       else
          !! In this case there is data in fileobj's root, so there is no need for the dummy data
          if(fileobj%is_root) allocate(global_buf_r8_kind(i_glob, j_glob))
          allocate(local_buf_r8_kind(size(vdata,1), size(vdata,2)))
          local_buf_r8_kind = vdata
       endif

       call mpp_gather(isc+i_add, iec+i_add, jsc+j_add, jec+j_add, bc_info%pelist, &
                          local_buf_r8_kind(i1:i2,j1:j2), &
                          global_buf_r8_kind, fileobj%is_root)
       deallocate(local_buf_r8_kind)
       !> If you are on fileobj's root, calculate the checksum and save the gathered data in a buffer
       if (fileobj%is_root) then
           chksum_val = mpp_chksum(global_buf_r8_kind(1:i_glob,1:j_glob), (/mpp_pe()/))
           allocate(bc_info%globaldata2d_r8(i_glob, j_glob))
           bc_info%globaldata2d_r8=global_buf_r8_kind(1:i_glob,1:j_glob)
           deallocate(global_buf_r8_kind)
       endif
     class default
         call error("gather_data_bc_2d: unsupported type. Currently only r8_kind and r4_kinds are supported")
   end select

   !> Save the checksum, so you can write it later
   if (fileobj%is_root) then
      chksum = ""
      write(chksum, "(Z16)") chksum_val
      bc_info%chksum = chksum
   endif

end subroutine gather_data_bc_2d

!> @brief gathers the 2d vdata from all of the relevant pes into the root_pe and saves it to a
!! buffer.
subroutine gather_data_bc_3d(fileobj, vdata, bc_info)
  class(FmsNetcdfFile_t), intent(inout) :: fileobj !< Fms2io netcdf fileobj
  class(*), dimension(:,:,:), intent(in) :: vdata !< Data to be gather
  type(bc_information), intent(inout) :: bc_info !< information about the boundary condition variable

  integer :: i_glob !< Size of the global domain in x
  integer :: j_glob !< Size of the global domain in y
  integer :: k_glob !< Size of the z axis

  integer :: isc, iec, jsc, jec !< current PE's indices
  integer :: i1, i2, j1, j2 !< current PE's indices relative to the global domain
  integer :: i_add, j_add !< current PE's shift

  real(kind=r4_kind), dimension(:,:,:), allocatable, target :: global_buf_r4_kind !< buffer with a data gathered
  real(kind=r4_kind), dimension(:,:,:), allocatable, target :: local_buf_r4_kind !< current PEs data
  real(kind=r8_kind), dimension(:,:,:), allocatable, target :: global_buf_r8_kind !< buffer with a data gathered
  real(kind=r8_kind), dimension(:,:,:), allocatable, target :: local_buf_r8_kind !< current PEs data

  integer(kind=i8_kind) :: chksum_val !< Checksum value calculated from the gathered data
  character(len=32) :: chksum !< Cheksum value converted to a string

  !> Set the indices
  isc = bc_info%indices(1)
  iec = bc_info%indices(2)
  jsc = bc_info%indices(3)
  jec = bc_info%indices(4)

  !> This is the section of the PE that will actually be added to the global_buffer
  i1 = 1 + bc_info%x_halo
  i2 = i1 + (iec-isc)
  j1 = 1 + bc_info%y_halo
  j2 = j1 + (jec-jsc)

  !> Set up index shifts for global array
  i_add = bc_info%ishift
  j_add = bc_info%jshift

  !> Allocate a global_buffer that will be written
  if (fileobj%is_root) then
        i_glob = bc_info%global_size(1)
        j_glob = bc_info%global_size(2)
  endif

  k_glob=bc_info%global_size(3)
  !> Gather the data and calculate the checksum for the resulting array.
  select type(vdata)
    type is (real(kind=r4_kind))
       !> If the fileobj's root does not have data for this variable
       if (fileobj%is_root .and. .not. bc_info%data_on_file_root) then
          !> Allocate global_buf_r8_kind to be one size bigger, global_buf_r8_kind(i_glob+1,,:) is just dummy data
          allocate(global_buf_r4_kind(i_glob+1, j_glob, bc_info%global_size(3)))
          !> Allocate a temp local buffer to the fileobj's root. This is needed because the data needs to be send
          !! to the fileobj's root, but because fileobj's root doesn't have any data, we just create dummy data and
          !! not write it later
          allocate(local_buf_r4_kind(1,1,1))
          local_buf_r4_kind = 0.
          isc = 1+i_glob; i_add=0; iec=1+i_glob; jsc=j_glob; j_add=0; jec=j_glob
          i1=1; i2=1; j1=1; j2=1
       else
          !! In this case there is data in fileobj's root, so there is no need for the dummy data
          if(fileobj%is_root) allocate(global_buf_r4_kind(i_glob, j_glob, k_glob))
          allocate(local_buf_r4_kind(size(vdata,1), size(vdata,2), size(vdata,3)))
          local_buf_r4_kind = vdata
       endif

       call mpp_gather(isc+i_add, iec+i_add, jsc+j_add, jec+j_add, k_glob, bc_info%pelist, &
                          local_buf_r4_kind(i1:i2,j1:j2,:), &
                          global_buf_r4_kind, fileobj%is_root)
       deallocate(local_buf_r4_kind)
       !> If you are on fileobj's root, calculate the checksum and save the gathered data in a buffer
       if (fileobj%is_root) then
           chksum_val = mpp_chksum(global_buf_r4_kind(1:i_glob,1:j_glob, :), (/mpp_pe()/))
           allocate(bc_info%globaldata3d_r4(i_glob, j_glob, bc_info%global_size(3)))
           bc_info%globaldata3d_r4=global_buf_r4_kind(1:i_glob,1:j_glob,:)
           deallocate(global_buf_r4_kind)
       endif
    type is (real(kind=r8_kind))
       !> If the fileobj's root does not have data for this variable
       if (fileobj%is_root .and. .not. bc_info%data_on_file_root) then
          !> Allocate global_buf_r8_kind to be one size bigger, global_buf_r8_kind(i_glob+1,,:) is just dummy data
          allocate(global_buf_r8_kind(i_glob+1, j_glob, bc_info%global_size(3)))
          !> Allocate a temp local buffer to the fileobj's root. This is needed because the data needs to be send
          !! to the fileobj's root, but because fileobj's root doesn't have any data, we just create dummy data and
          !! not write it later
          allocate(local_buf_r8_kind(1,1,1))
          local_buf_r8_kind = 0.
          isc = 1+i_glob; i_add=0; iec=1+i_glob; jsc=j_glob; j_add=0; jec=j_glob
          i1=1; i2=1; j1=1; j2=1
       else
          !! In this case there is data in fileobj's root, so there is no need for the dummy data
          if(fileobj%is_root) allocate(global_buf_r8_kind(i_glob, j_glob, k_glob))
          allocate(local_buf_r8_kind(size(vdata,1), size(vdata,2), size(vdata,3)))
          local_buf_r8_kind = vdata
       endif

       call mpp_gather(isc+i_add, iec+i_add, jsc+j_add, jec+j_add, k_glob, bc_info%pelist, &
                          local_buf_r8_kind(i1:i2,j1:j2,:), &
                          global_buf_r8_kind, fileobj%is_root)
       deallocate(local_buf_r8_kind)
       !> If you are on fileobj's root, calculate the checksum and save the gathered data in a buffer
       if (fileobj%is_root) then
           chksum_val = mpp_chksum(global_buf_r8_kind(1:i_glob,1:j_glob, :), (/mpp_pe()/))
           allocate(bc_info%globaldata3d_r8(i_glob, j_glob, bc_info%global_size(3)))
           bc_info%globaldata3d_r8=global_buf_r8_kind(1:i_glob,1:j_glob,:)
           deallocate(global_buf_r8_kind)
       endif
     class default
         call error("gather_data_bc_3d: unsupported type. Currently only r8_kind and r4_kinds are supported")
   end select

   !> Save the checksum, so you can write it later
   if (fileobj%is_root) then
      chksum = ""
      write(chksum, "(Z16)") chksum_val
      bc_info%chksum = chksum
   endif

end subroutine gather_data_bc_3d
