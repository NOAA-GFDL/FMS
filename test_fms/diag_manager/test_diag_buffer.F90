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

!> @brief This program tests the output buffer functionality
program test_diag_buffer
#ifdef use_yaml

    use fms_diag_output_buffer_mod, only: fmsDiagOutputBuffer_type
    use platform_mod,               only: r8_kind, r4_kind, i8_kind, i4_kind
    use fms_mod,                    only: string, fms_init, fms_end
    use mpp_mod,                    only: mpp_error, FATAL
    use diag_data_mod,              only: i4, i8, r4, r8, time_none, EMPTY

    implicit none

    type(fmsDiagOutputBuffer_type) :: buffobj(6)       !< Dummy output buffers
    integer                   :: buff_sizes(5)    !< Size of the buffer for each dimension
    class(*),allocatable      :: p_val(:,:,:,:,:) !< Dummy variable to get the data
    integer                   :: i, j             !< For do loops
    real(r8_kind)             :: r8_data          !< Dummy r8 data
    real(r4_kind)             :: r4_data          !< Dummy r4 data
    integer(i8_kind)          :: i8_data          !< Dummy i8 data
    integer(i4_kind)          :: i4_data          !< Dummy i4 data
    character(len=4)          :: fname = 'test'   !< Dummy name for error messages

    call fms_init

    !< Test the r8_buffer
    buff_sizes = 1
    do i=0, 4
        buff_sizes(i+1) = i+5
        call buffobj(i+1)%allocate_buffer(r8_data, i, buff_sizes, .false., fname, 1)
        call buffobj(i+1)%initialize_buffer(time_none, fname)
        call buffobj(i+1)%get_buffer(p_val, fname)
        select type(p_val)
        type is (real(kind=r8_kind))
          if (any(p_val .ne. real(EMPTY, kind=r8_kind))) &
            call mpp_error(FATAL, "r8_buffer:: The "//string(i)//"d buffer was not initialized to the correct value")
        do j = 1, 4
          if (size(p_val, j) .ne. buff_sizes(j)) then
            call mpp_error(FATAL, "r8_buffer:: The "//string(i)//"d buffer was not allocated to the correct size")
          endif
        enddo
        class default
          call mpp_error(FATAL, "r8_buffer:: The "//string(i)//"d buffer was not allocated to the correct type")
        end select
        deallocate(p_val)
        call buffobj(i+1)%flush_buffer()
    end do

    !< Test the r4_buffer
    buff_sizes = 1
    do i=0, 4
        buff_sizes(i+1) = i+5
        call buffobj(i+1)%allocate_buffer(r4_data, i, buff_sizes, .false., fname, 1)
        call buffobj(i+1)%initialize_buffer(time_none, fname)
        call buffobj(i+1)%get_buffer(p_val, fname)
        select type(p_val)
        type is (real(kind=r4_kind))
          if (any(p_val .ne. real(EMPTY, kind=r4_kind))) &
            call mpp_error(FATAL, "r4_buffer:: The "//string(i)//"d buffer was not initialized to the correct value")
        do j = 1, 4
          if (size(p_val, j) .ne. buff_sizes(j)) &
            call mpp_error(FATAL, "r4_buffer:: The "//string(i)//"d buffer was not allocated to the correct size")
        enddo
        class default
          call mpp_error(FATAL, "r4_buffer:: The "//string(i)//"d buffer was not allocated to the correct type")
        end select
        deallocate(p_val)
        call buffobj(i+1)%flush_buffer()
    end do

    !< Test the i8_buffer
    buff_sizes = 1
    do i=0, 4
        buff_sizes(i+1) = i+5
        call buffobj(i+1)%allocate_buffer(i8_data, i, buff_sizes, .false., fname, 1)
        call buffobj(i+1)%initialize_buffer(time_none, fname)
        call buffobj(i+1)%get_buffer(p_val, fname)
        select type(p_val)
        type is (integer(kind=i8_kind))
          if (any(p_val .ne. int(EMPTY, kind=i8_kind))) &
            call mpp_error(FATAL, "i8_buffer:: The "//string(i)//"d buffer was not initialized to the correct value")
        do j = 1, 4
          if (size(p_val, j) .ne. buff_sizes(j)) &
            call mpp_error(FATAL, "i8_buffer:: The "//string(i)//"d buffer was not allocated to the correct size")
        enddo
        class default
          call mpp_error(FATAL, "i8_buffer:: The "//string(i)//"d buffer was not allocated to the correct type")
        end select
        deallocate(p_val)
        call buffobj(i+1)%flush_buffer()
    end do

    !< Test the i4_buffer
    buff_sizes = 1
    do i=0, 4
        buff_sizes(i+1) = i+5
        call buffobj(i+1)%allocate_buffer(i4_data, i, buff_sizes, .false., fname, 1)
        call buffobj(i+1)%initialize_buffer(time_none, fname)
        call buffobj(i+1)%get_buffer(p_val, fname)
        select type(p_val)
        type is (integer(kind=i4_kind))
          if (any(p_val .ne. int(EMPTY, kind=i4_kind))) &
            call mpp_error(FATAL, "i4_buffer:: The "//string(i)//"d buffer was not initialized to the correct value")
        do j = 1, 4
          if (size(p_val, j) .ne. buff_sizes(j)) &
            call mpp_error(FATAL, "i4_buffer:: The "//string(i)//"d buffer was not allocated to the correct size")
        enddo
        class default
          call mpp_error(FATAL, "i4_buffer:: The "//string(i)//"d buffer was not allocated to the correct type")
        end select
        deallocate(p_val)
        call buffobj(i+1)%flush_buffer()
    end do

    call fms_end()
#endif
end program
