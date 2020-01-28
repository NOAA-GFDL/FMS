!> @brief Allocate character arrays using an input array of sizes.
subroutine allocate_array_char_1d(buf, sizes)

  character(len=*), dimension(:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(1), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1)))
end subroutine allocate_array_char_1d


!> @brief Allocate character arrays using an input array of sizes.
subroutine allocate_array_char_2d(buf, sizes)

  character(len=*), dimension(:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(2), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1), sizes(2)))
end subroutine allocate_array_char_2d


!> @brief Allocate character arrays using an input array of sizes.
subroutine allocate_array_char_3d(buf, sizes)

  character(len=*), dimension(:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(3), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1), sizes(2), sizes(3)))
end subroutine allocate_array_char_3d


!> @brief Allocate character arrays using an input array of sizes.
subroutine allocate_array_char_4d(buf, sizes)

  character(len=*), dimension(:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(4), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1), sizes(2), sizes(3), sizes(4)))
end subroutine allocate_array_char_4d


!> @brief Allocate character arrays using an input array of sizes.
subroutine allocate_array_char_5d(buf, sizes)

  character(len=*), dimension(:,:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(5), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1), sizes(2), sizes(3), sizes(4), sizes(5)))
end subroutine allocate_array_char_5d


!> @brief Allocate character arrays using an input array of sizes.
subroutine allocate_array_char_6d(buf, sizes)

  character(len=*), dimension(:,:,:,:,:,:), allocatable, intent(inout) :: buf !< Array that will be allocated.
  integer, dimension(6), intent(in) :: sizes !< Array of dimension sizes.

  if (allocated(buf)) then
    deallocate(buf)
  endif
  allocate(buf(sizes(1), sizes(2), sizes(3), sizes(4), sizes(5), sizes(6)))
end subroutine allocate_array_char_6d


