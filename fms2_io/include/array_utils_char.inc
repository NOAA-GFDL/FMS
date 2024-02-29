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
!> @brief Character array routines used in @ref fms_io_utils_mod

!> @addtogroup fms_io_utils_mod
!> @{

!> @brief Allocate character arrays using an input array of sizes.
subroutine allocate_array_char_1d(buf, sizes, initialize)

  character(len=*), dimension(:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(1), intent(in) :: sizes !< Array of dimension sizes.
  logical, intent(in), optional :: initialize !< Optional argument when true will initialize with a blank string.

  logical :: init !< local variable for initialize
  integer :: i, c !< for looping

  init = .false.
  if (present(initialize)) init = initialize

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1)))

  if (init) then
    do i = 1, sizes(1)
      do c = 1, len(buf(i))
        buf(i)(c:c) = " "
      enddo
    enddo
  endif

end subroutine allocate_array_char_1d


!> @brief Allocate character arrays using an input array of sizes.
subroutine allocate_array_char_2d(buf, sizes, initialize)

  character(len=*), dimension(:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(2), intent(in) :: sizes !< Array of dimension sizes.
  logical, intent(in), optional :: initialize !< Optional argument when true will initialize with a blank string.

  logical :: init !< local variable for initialize
  integer :: i, j, c !< for looping

  init = .false.
  if (present(initialize)) init = initialize

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1), sizes(2)))

  if (init) then
    do j = 1, sizes(2)
      do i = 1, sizes(1)
        do c = 1, len(buf(i,j))
          buf(i,j)(c:c) = " "
        enddo
      enddo
    enddo
  endif

end subroutine allocate_array_char_2d


!> @brief Allocate character arrays using an input array of sizes.
subroutine allocate_array_char_3d(buf, sizes, initialize)

  character(len=*), dimension(:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(3), intent(in) :: sizes !< Array of dimension sizes.
  logical, intent(in), optional :: initialize !< Optional argument when true will initialize with a blank string.

  logical :: init !< local variable for initialize
  integer :: i, j, k, c !< for looping

  init = .false.
  if (present(initialize)) init = initialize

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1), sizes(2), sizes(3)))

  if (init) then
    do k = 1, sizes(3)
      do j = 1, sizes(2)
        do i = 1, sizes(1)
          do c = 1, len(buf(i,j,k))
            buf(i,j,k)(c:c) = " "
          enddo
        enddo
      enddo
    enddo
  endif

end subroutine allocate_array_char_3d


!> @brief Allocate character arrays using an input array of sizes.
subroutine allocate_array_char_4d(buf, sizes, initialize)

  character(len=*), dimension(:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(4), intent(in) :: sizes !< Array of dimension sizes.
  logical, intent(in), optional :: initialize !< Optional argument when true will initialize with a blank string.

  logical :: init !< local variable for initialize
  integer :: i, j, k, l, c !< for looping

  init = .false.
  if (present(initialize)) init = initialize

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1), sizes(2), sizes(3), sizes(4)))

  if (init) then
    do l = 1, sizes(4)
      do k = 1, sizes(3)
        do j = 1, sizes(2)
          do i = 1, sizes(1)
            do c = 1, len(buf(i,j,k,l))
              buf(i,j,k,l)(c:c) = " "
            enddo
          enddo
        enddo
      enddo
    enddo
  endif
end subroutine allocate_array_char_4d


!> @brief Allocate character arrays using an input array of sizes.
subroutine allocate_array_char_5d(buf, sizes, initialize)

  character(len=*), dimension(:,:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(5), intent(in) :: sizes !< Array of dimension sizes.
  logical, intent(in), optional :: initialize !< Optional argument when true will initialize with a blank string.

  logical :: init !< local variable for initialize
  integer :: i, j, k, l, m, c !< for looping

  init = .false.
  if (present(initialize)) init = initialize

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1), sizes(2), sizes(3), sizes(4), sizes(5)))

  if (init) then
    do m = 1, sizes(5)
      do l = 1, sizes(4)
        do k = 1, sizes(3)
          do j = 1, sizes(2)
            do i = 1, sizes(1)
              do c = 1, len(buf(i,j,k,l,m))
                buf(i,j,k,l,m)(c:c) = " "
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  endif
end subroutine allocate_array_char_5d


!> @brief Allocate character arrays using an input array of sizes.
subroutine allocate_array_char_6d(buf, sizes, initialize)

  character(len=*), dimension(:,:,:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(6), intent(in) :: sizes !< Array of dimension sizes.
  logical, intent(in), optional :: initialize !< Optional argument when true will initialize with a blank string.

  logical :: init !< local variable for initialize
  integer :: i, j, k, l, m, n, c !< for looping

  init = .false.
  if (present(initialize)) init = initialize

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1), sizes(2), sizes(3), sizes(4), sizes(5), sizes(6)))

  if (init) then
    do n = 1, sizes(6)
      do m = 1, sizes(5)
        do l = 1, sizes(4)
          do k = 1, sizes(3)
            do j = 1, sizes(2)
              do i = 1, sizes(1)
                do c = 1, len(buf(i,j,k,l,m,n))
                  buf(i,j,k,l,m,n)(c:c) = " "
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  endif
end subroutine allocate_array_char_6d
!> @}
