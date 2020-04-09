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

!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_int32_1d(buf, sizes)

  integer(kind=int32), dimension(:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(1), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1)))
end subroutine allocate_array_int32_1d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_int32_1d(section, array, s, c)

  integer(kind=int32), dimension(:), intent(in) :: section !< Section to be inserted.
  integer(kind=int32), dimension(:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(1), intent(in) :: s !< Array of starting indices.
  integer, dimension(1), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ) = section(:)
end subroutine put_array_section_int32_1d


!> @brief Get a section of larger array.
subroutine get_array_section_int32_1d(section, array, s, c)

  integer(kind=int32), dimension(:), intent(inout) :: section !< Section to be extracted.
  integer(kind=int32), dimension(:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(1), intent(in) :: s !< Array of starting indices.
  integer, dimension(1), intent(in) :: c !< Array of sizes.

  section(:) = array(s(1):s(1)+c(1)-1 )
end subroutine get_array_section_int32_1d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_int32_2d(buf, sizes)

  integer(kind=int32), dimension(:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(2), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2)))
end subroutine allocate_array_int32_2d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_int32_2d(section, array, s, c)

  integer(kind=int32), dimension(:,:), intent(in) :: section !< Section to be inserted.
  integer(kind=int32), dimension(:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(2), intent(in) :: s !< Array of starting indices.
  integer, dimension(2), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ) = section(:,:)
end subroutine put_array_section_int32_2d


!> @brief Get a section of larger array.
subroutine get_array_section_int32_2d(section, array, s, c)

  integer(kind=int32), dimension(:,:), intent(inout) :: section !< Section to be extracted.
  integer(kind=int32), dimension(:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(2), intent(in) :: s !< Array of starting indices.
  integer, dimension(2), intent(in) :: c !< Array of sizes.

  section(:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 )
end subroutine get_array_section_int32_2d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_int32_3d(buf, sizes)

  integer(kind=int32), dimension(:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(3), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2),sizes(3)))
end subroutine allocate_array_int32_3d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_int32_3d(section, array, s, c)

  integer(kind=int32), dimension(:,:,:), intent(in) :: section !< Section to be inserted.
  integer(kind=int32), dimension(:,:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(3), intent(in) :: s !< Array of starting indices.
  integer, dimension(3), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ) = section(:,:,:)
end subroutine put_array_section_int32_3d


!> @brief Get a section of larger array.
subroutine get_array_section_int32_3d(section, array, s, c)

  integer(kind=int32), dimension(:,:,:), intent(inout) :: section !< Section to be extracted.
  integer(kind=int32), dimension(:,:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(3), intent(in) :: s !< Array of starting indices.
  integer, dimension(3), intent(in) :: c !< Array of sizes.

  section(:,:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 )
end subroutine get_array_section_int32_3d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_int32_4d(buf, sizes)

  integer(kind=int32), dimension(:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(4), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2),sizes(3),sizes(4)))
end subroutine allocate_array_int32_4d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_int32_4d(section, array, s, c)

  integer(kind=int32), dimension(:,:,:,:), intent(in) :: section !< Section to be inserted.
  integer(kind=int32), dimension(:,:,:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(4), intent(in) :: s !< Array of starting indices.
  integer, dimension(4), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 ) = section(:,:,:,:)
end subroutine put_array_section_int32_4d


!> @brief Get a section of larger array.
subroutine get_array_section_int32_4d(section, array, s, c)

  integer(kind=int32), dimension(:,:,:,:), intent(inout) :: section !< Section to be extracted.
  integer(kind=int32), dimension(:,:,:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(4), intent(in) :: s !< Array of starting indices.
  integer, dimension(4), intent(in) :: c !< Array of sizes.

  section(:,:,:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 )
end subroutine get_array_section_int32_4d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_int32_5d(buf, sizes)

  integer(kind=int32), dimension(:,:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(5), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2),sizes(3),sizes(4),sizes(5)))
end subroutine allocate_array_int32_5d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_int32_5d(section, array, s, c)

  integer(kind=int32), dimension(:,:,:,:,:), intent(in) :: section !< Section to be inserted.
  integer(kind=int32), dimension(:,:,:,:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(5), intent(in) :: s !< Array of starting indices.
  integer, dimension(5), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 ,s(5):s(5)+c(5)-1 ) = section(:,:,:,:,:)
end subroutine put_array_section_int32_5d


!> @brief Get a section of larger array.
subroutine get_array_section_int32_5d(section, array, s, c)

  integer(kind=int32), dimension(:,:,:,:,:), intent(inout) :: section !< Section to be extracted.
  integer(kind=int32), dimension(:,:,:,:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(5), intent(in) :: s !< Array of starting indices.
  integer, dimension(5), intent(in) :: c !< Array of sizes.

  section(:,:,:,:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 ,s(5):s(5)+c(5)-1 )
end subroutine get_array_section_int32_5d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_int64_1d(buf, sizes)

  integer(kind=int64), dimension(:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(1), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1)))
end subroutine allocate_array_int64_1d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_int64_1d(section, array, s, c)

  integer(kind=int64), dimension(:), intent(in) :: section !< Section to be inserted.
  integer(kind=int64), dimension(:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(1), intent(in) :: s !< Array of starting indices.
  integer, dimension(1), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ) = section(:)
end subroutine put_array_section_int64_1d


!> @brief Get a section of larger array.
subroutine get_array_section_int64_1d(section, array, s, c)

  integer(kind=int64), dimension(:), intent(inout) :: section !< Section to be extracted.
  integer(kind=int64), dimension(:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(1), intent(in) :: s !< Array of starting indices.
  integer, dimension(1), intent(in) :: c !< Array of sizes.

  section(:) = array(s(1):s(1)+c(1)-1 )
end subroutine get_array_section_int64_1d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_int64_2d(buf, sizes)

  integer(kind=int64), dimension(:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(2), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2)))
end subroutine allocate_array_int64_2d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_int64_2d(section, array, s, c)

  integer(kind=int64), dimension(:,:), intent(in) :: section !< Section to be inserted.
  integer(kind=int64), dimension(:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(2), intent(in) :: s !< Array of starting indices.
  integer, dimension(2), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ) = section(:,:)
end subroutine put_array_section_int64_2d


!> @brief Get a section of larger array.
subroutine get_array_section_int64_2d(section, array, s, c)

  integer(kind=int64), dimension(:,:), intent(inout) :: section !< Section to be extracted.
  integer(kind=int64), dimension(:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(2), intent(in) :: s !< Array of starting indices.
  integer, dimension(2), intent(in) :: c !< Array of sizes.

  section(:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 )
end subroutine get_array_section_int64_2d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_int64_3d(buf, sizes)

  integer(kind=int64), dimension(:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(3), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2),sizes(3)))
end subroutine allocate_array_int64_3d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_int64_3d(section, array, s, c)

  integer(kind=int64), dimension(:,:,:), intent(in) :: section !< Section to be inserted.
  integer(kind=int64), dimension(:,:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(3), intent(in) :: s !< Array of starting indices.
  integer, dimension(3), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ) = section(:,:,:)
end subroutine put_array_section_int64_3d


!> @brief Get a section of larger array.
subroutine get_array_section_int64_3d(section, array, s, c)

  integer(kind=int64), dimension(:,:,:), intent(inout) :: section !< Section to be extracted.
  integer(kind=int64), dimension(:,:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(3), intent(in) :: s !< Array of starting indices.
  integer, dimension(3), intent(in) :: c !< Array of sizes.

  section(:,:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 )
end subroutine get_array_section_int64_3d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_int64_4d(buf, sizes)

  integer(kind=int64), dimension(:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(4), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2),sizes(3),sizes(4)))
end subroutine allocate_array_int64_4d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_int64_4d(section, array, s, c)

  integer(kind=int64), dimension(:,:,:,:), intent(in) :: section !< Section to be inserted.
  integer(kind=int64), dimension(:,:,:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(4), intent(in) :: s !< Array of starting indices.
  integer, dimension(4), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 ) = section(:,:,:,:)
end subroutine put_array_section_int64_4d


!> @brief Get a section of larger array.
subroutine get_array_section_int64_4d(section, array, s, c)

  integer(kind=int64), dimension(:,:,:,:), intent(inout) :: section !< Section to be extracted.
  integer(kind=int64), dimension(:,:,:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(4), intent(in) :: s !< Array of starting indices.
  integer, dimension(4), intent(in) :: c !< Array of sizes.

  section(:,:,:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 )
end subroutine get_array_section_int64_4d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_int64_5d(buf, sizes)

  integer(kind=int64), dimension(:,:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(5), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2),sizes(3),sizes(4),sizes(5)))
end subroutine allocate_array_int64_5d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_int64_5d(section, array, s, c)

  integer(kind=int64), dimension(:,:,:,:,:), intent(in) :: section !< Section to be inserted.
  integer(kind=int64), dimension(:,:,:,:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(5), intent(in) :: s !< Array of starting indices.
  integer, dimension(5), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 ,s(5):s(5)+c(5)-1 ) = section(:,:,:,:,:)
end subroutine put_array_section_int64_5d


!> @brief Get a section of larger array.
subroutine get_array_section_int64_5d(section, array, s, c)

  integer(kind=int64), dimension(:,:,:,:,:), intent(inout) :: section !< Section to be extracted.
  integer(kind=int64), dimension(:,:,:,:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(5), intent(in) :: s !< Array of starting indices.
  integer, dimension(5), intent(in) :: c !< Array of sizes.

  section(:,:,:,:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 ,s(5):s(5)+c(5)-1 )
end subroutine get_array_section_int64_5d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_real32_1d(buf, sizes)

  real(kind=real32), dimension(:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(1), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1)))
end subroutine allocate_array_real32_1d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_real32_1d(section, array, s, c)

  real(kind=real32), dimension(:), intent(in) :: section !< Section to be inserted.
  real(kind=real32), dimension(:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(1), intent(in) :: s !< Array of starting indices.
  integer, dimension(1), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ) = section(:)
end subroutine put_array_section_real32_1d


!> @brief Get a section of larger array.
subroutine get_array_section_real32_1d(section, array, s, c)

  real(kind=real32), dimension(:), intent(inout) :: section !< Section to be extracted.
  real(kind=real32), dimension(:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(1), intent(in) :: s !< Array of starting indices.
  integer, dimension(1), intent(in) :: c !< Array of sizes.

  section(:) = array(s(1):s(1)+c(1)-1 )
end subroutine get_array_section_real32_1d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_real32_2d(buf, sizes)

  real(kind=real32), dimension(:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(2), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2)))
end subroutine allocate_array_real32_2d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_real32_2d(section, array, s, c)

  real(kind=real32), dimension(:,:), intent(in) :: section !< Section to be inserted.
  real(kind=real32), dimension(:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(2), intent(in) :: s !< Array of starting indices.
  integer, dimension(2), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ) = section(:,:)
end subroutine put_array_section_real32_2d


!> @brief Get a section of larger array.
subroutine get_array_section_real32_2d(section, array, s, c)

  real(kind=real32), dimension(:,:), intent(inout) :: section !< Section to be extracted.
  real(kind=real32), dimension(:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(2), intent(in) :: s !< Array of starting indices.
  integer, dimension(2), intent(in) :: c !< Array of sizes.

  section(:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 )
end subroutine get_array_section_real32_2d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_real32_3d(buf, sizes)

  real(kind=real32), dimension(:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(3), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2),sizes(3)))
end subroutine allocate_array_real32_3d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_real32_3d(section, array, s, c)

  real(kind=real32), dimension(:,:,:), intent(in) :: section !< Section to be inserted.
  real(kind=real32), dimension(:,:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(3), intent(in) :: s !< Array of starting indices.
  integer, dimension(3), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ) = section(:,:,:)
end subroutine put_array_section_real32_3d


!> @brief Get a section of larger array.
subroutine get_array_section_real32_3d(section, array, s, c)

  real(kind=real32), dimension(:,:,:), intent(inout) :: section !< Section to be extracted.
  real(kind=real32), dimension(:,:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(3), intent(in) :: s !< Array of starting indices.
  integer, dimension(3), intent(in) :: c !< Array of sizes.

  section(:,:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 )
end subroutine get_array_section_real32_3d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_real32_4d(buf, sizes)

  real(kind=real32), dimension(:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(4), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2),sizes(3),sizes(4)))
end subroutine allocate_array_real32_4d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_real32_4d(section, array, s, c)

  real(kind=real32), dimension(:,:,:,:), intent(in) :: section !< Section to be inserted.
  real(kind=real32), dimension(:,:,:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(4), intent(in) :: s !< Array of starting indices.
  integer, dimension(4), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 ) = section(:,:,:,:)
end subroutine put_array_section_real32_4d


!> @brief Get a section of larger array.
subroutine get_array_section_real32_4d(section, array, s, c)

  real(kind=real32), dimension(:,:,:,:), intent(inout) :: section !< Section to be extracted.
  real(kind=real32), dimension(:,:,:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(4), intent(in) :: s !< Array of starting indices.
  integer, dimension(4), intent(in) :: c !< Array of sizes.

  section(:,:,:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 )
end subroutine get_array_section_real32_4d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_real32_5d(buf, sizes)

  real(kind=real32), dimension(:,:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(5), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2),sizes(3),sizes(4),sizes(5)))
end subroutine allocate_array_real32_5d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_real32_5d(section, array, s, c)

  real(kind=real32), dimension(:,:,:,:,:), intent(in) :: section !< Section to be inserted.
  real(kind=real32), dimension(:,:,:,:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(5), intent(in) :: s !< Array of starting indices.
  integer, dimension(5), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 ,s(5):s(5)+c(5)-1 ) = section(:,:,:,:,:)
end subroutine put_array_section_real32_5d


!> @brief Get a section of larger array.
subroutine get_array_section_real32_5d(section, array, s, c)

  real(kind=real32), dimension(:,:,:,:,:), intent(inout) :: section !< Section to be extracted.
  real(kind=real32), dimension(:,:,:,:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(5), intent(in) :: s !< Array of starting indices.
  integer, dimension(5), intent(in) :: c !< Array of sizes.

  section(:,:,:,:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 ,s(5):s(5)+c(5)-1 )
end subroutine get_array_section_real32_5d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_real64_1d(buf, sizes)

  real(kind=real64), dimension(:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(1), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1)))
end subroutine allocate_array_real64_1d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_real64_1d(section, array, s, c)

  real(kind=real64), dimension(:), intent(in) :: section !< Section to be inserted.
  real(kind=real64), dimension(:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(1), intent(in) :: s !< Array of starting indices.
  integer, dimension(1), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ) = section(:)
end subroutine put_array_section_real64_1d


!> @brief Get a section of larger array.
subroutine get_array_section_real64_1d(section, array, s, c)

  real(kind=real64), dimension(:), intent(inout) :: section !< Section to be extracted.
  real(kind=real64), dimension(:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(1), intent(in) :: s !< Array of starting indices.
  integer, dimension(1), intent(in) :: c !< Array of sizes.

  section(:) = array(s(1):s(1)+c(1)-1 )
end subroutine get_array_section_real64_1d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_real64_2d(buf, sizes)

  real(kind=real64), dimension(:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(2), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2)))
end subroutine allocate_array_real64_2d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_real64_2d(section, array, s, c)

  real(kind=real64), dimension(:,:), intent(in) :: section !< Section to be inserted.
  real(kind=real64), dimension(:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(2), intent(in) :: s !< Array of starting indices.
  integer, dimension(2), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ) = section(:,:)
end subroutine put_array_section_real64_2d


!> @brief Get a section of larger array.
subroutine get_array_section_real64_2d(section, array, s, c)

  real(kind=real64), dimension(:,:), intent(inout) :: section !< Section to be extracted.
  real(kind=real64), dimension(:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(2), intent(in) :: s !< Array of starting indices.
  integer, dimension(2), intent(in) :: c !< Array of sizes.

  section(:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 )
end subroutine get_array_section_real64_2d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_real64_3d(buf, sizes)

  real(kind=real64), dimension(:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(3), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2),sizes(3)))
end subroutine allocate_array_real64_3d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_real64_3d(section, array, s, c)

  real(kind=real64), dimension(:,:,:), intent(in) :: section !< Section to be inserted.
  real(kind=real64), dimension(:,:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(3), intent(in) :: s !< Array of starting indices.
  integer, dimension(3), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ) = section(:,:,:)
end subroutine put_array_section_real64_3d


!> @brief Get a section of larger array.
subroutine get_array_section_real64_3d(section, array, s, c)

  real(kind=real64), dimension(:,:,:), intent(inout) :: section !< Section to be extracted.
  real(kind=real64), dimension(:,:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(3), intent(in) :: s !< Array of starting indices.
  integer, dimension(3), intent(in) :: c !< Array of sizes.

  section(:,:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 )
end subroutine get_array_section_real64_3d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_real64_4d(buf, sizes)

  real(kind=real64), dimension(:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(4), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2),sizes(3),sizes(4)))
end subroutine allocate_array_real64_4d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_real64_4d(section, array, s, c)

  real(kind=real64), dimension(:,:,:,:), intent(in) :: section !< Section to be inserted.
  real(kind=real64), dimension(:,:,:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(4), intent(in) :: s !< Array of starting indices.
  integer, dimension(4), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 ) = section(:,:,:,:)
end subroutine put_array_section_real64_4d


!> @brief Get a section of larger array.
subroutine get_array_section_real64_4d(section, array, s, c)

  real(kind=real64), dimension(:,:,:,:), intent(inout) :: section !< Section to be extracted.
  real(kind=real64), dimension(:,:,:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(4), intent(in) :: s !< Array of starting indices.
  integer, dimension(4), intent(in) :: c !< Array of sizes.

  section(:,:,:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 )
end subroutine get_array_section_real64_4d


!> @brief Allocate arrays using an input array of sizes.
subroutine allocate_array_real64_5d(buf, sizes)

  real(kind=real64), dimension(:,:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(5), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1),sizes(2),sizes(3),sizes(4),sizes(5)))
end subroutine allocate_array_real64_5d


!> @brief Put a section of an array into a larger array.
subroutine put_array_section_real64_5d(section, array, s, c)

  real(kind=real64), dimension(:,:,:,:,:), intent(in) :: section !< Section to be inserted.
  real(kind=real64), dimension(:,:,:,:,:), intent(inout) :: array !< Array to insert the section in.
  integer, dimension(5), intent(in) :: s !< Array of starting indices.
  integer, dimension(5), intent(in) :: c !< Array of sizes.

  array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 ,s(5):s(5)+c(5)-1 ) = section(:,:,:,:,:)
end subroutine put_array_section_real64_5d


!> @brief Get a section of larger array.
subroutine get_array_section_real64_5d(section, array, s, c)

  real(kind=real64), dimension(:,:,:,:,:), intent(inout) :: section !< Section to be extracted.
  real(kind=real64), dimension(:,:,:,:,:), intent(in) :: array !< Array to extract the section from.
  integer, dimension(5), intent(in) :: s !< Array of starting indices.
  integer, dimension(5), intent(in) :: c !< Array of sizes.

  section(:,:,:,:,:) = array(s(1):s(1)+c(1)-1 ,s(2):s(2)+c(2)-1 ,s(3):s(3)+c(3)-1 ,s(4):s(4)+c(4)-1 ,s(5):s(5)+c(5)-1 )
end subroutine get_array_section_real64_5d


