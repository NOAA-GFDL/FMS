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

!> @brief Utilities used in multiple test
module testing_utils
  use platform_mod,      only: r8_kind
  use mpp_mod,           only: mpp_error, FATAL
  private

  public :: allocate_buffer, permute, check_perm
  public :: test_normal, test_openmp, test_halos
  public :: no_mask, logical_mask, real_mask

  integer, parameter :: test_normal = 0  !< sending a buffer in the compute domain
  integer, parameter :: test_openmp = 1  !< sending a buffer in the compute domain but with blocking
  integer, parameter :: test_halos = 2   !< sending a buffer in the data domain (i.e with halos)
  integer, parameter :: no_mask = 0      !< Not using a mask
  integer, parameter :: logical_mask = 1 !< Using a logical mask
  integer, parameter :: real_mask = 2    !< Using a real mask

  interface permute
    module procedure permute_2d
    module procedure permute_3d
  end interface permute

  interface check_perm
    module procedure check_perm_2d
    module procedure check_perm_3d
  end interface check_perm

  contains

  !> @brief Allocate the output buffer based on the starting/ending indices
  !! @return output buffer set to -999_r8_kind
  function allocate_buffer(is, ie, js, je, k, l) &
    result(buffer)
    integer, intent(in) :: is !< Starting x index
    integer, intent(in) :: ie !< Ending x index
    integer, intent(in) :: js !< Starting y index
    integer, intent(in) :: je !< Ending y index
    integer, intent(in) :: k  !< Number of points in the 4th dimension
    integer, intent(in) :: l  !< Number of points in the 5th dimension
    real(kind=r8_kind), allocatable :: buffer(:,:,:,:)

    allocate(buffer(is:ie, js:je, 1:k, 1:l))
    buffer = -999_r8_kind
  end function allocate_buffer

  !> @brief Apply a 2D axis permutation to an array.
  !! @return permuted array with axes reordered according to order
  function permute_2d(src, order) result(dst)
    real(r8_kind), intent(in)  :: src(:,:)
    integer,       intent(in)  :: order(2)
    real(r8_kind), allocatable :: dst(:,:)

    integer :: i, j
    integer :: idx(2)

    allocate(dst(size(src,order(1)), size(src,order(2))))

    do j = 1, size(src,2)
       do i = 1, size(src,1)
         idx = [i,j]
         dst(idx(order(1)), idx(order(2))) = src(i,j)
       end do
    end do
  end function permute_2d

  !> @brief Apply a 3D axis permutation to an array.
  !! @return permuted array with axes reordered according to order
  function permute_3d(src, order) result(dst)
    real(r8_kind), intent(in)  :: src(:,:,:)
    integer,       intent(in)  :: order(3)
    real(r8_kind), allocatable :: dst(:,:,:)

    integer :: i, j, k
    integer :: idx(3)

    allocate(dst( size(src,order(1)), size(src,order(2)), size(src,order(3)) ))

    do k = 1, size(src,3)
       do j = 1, size(src,2)
          do i = 1, size(src,1)
             idx = [i,j,k]
             dst(idx(order(1)), idx(order(2)), idx(order(3))) = src(i,j,k)
          enddo
       enddo
    enddo
  end function permute_3d

  !> @brief Verify correctness of a 2D axis permutation.
  !! Aborts if shape or values do not match expected permutation.
  subroutine check_perm_2d(var, var_perm, order)
    real(kind=r8_kind), intent(in) :: var(:,:)      ! canonical (x,y)
    real(kind=r8_kind), intent(in) :: var_perm(:,:) ! permuted
    integer,            intent(in) :: order(2)

    integer :: i, j
    integer :: idx(2)

    ! Check shape consistency
    if ( size(var,order(1)) /= size(var_perm,1) .or. &
         size(var,order(2)) /= size(var_perm,2) ) then
       call mpp_error(FATAL, "check_perm_2d: dimension mismatch")
    endif

    do j = 1, size(var,2)
       do i = 1, size(var,1)
          idx = [i,j]

          if (abs(var(i,j) - var_perm(idx(order(1)), idx(order(2)))) > 0) then
             print *, "perm mismatch at (x,y)=", i, j, "order=", order, &
                         " var  =", var(i,j), &
                         " perm =", var_perm(idx(order(1)), idx(order(2)))
             call mpp_error(FATAL, "check_perm_2d failed")
          endif

       end do
    end do
  end subroutine check_perm_2d

  !> @brief Verify correctness of a 3D axis permutation.
  !! Aborts if shape or values do not match expected permutation.
  subroutine check_perm_3d(var, var_perm, order)
    real(kind=r8_kind), intent(in) :: var(:,:,:)      ! canonical (x,y,z)
    real(kind=r8_kind), intent(in) :: var_perm(:,:,:) ! permuted
    integer,            intent(in) :: order(3)

    integer :: i, j, k
    integer :: idx(3)

    ! Check shape consistency
    if ( size(var,order(1)) /= size(var_perm,1) .or. &
         size(var,order(2)) /= size(var_perm,2) .or. &
         size(var,order(3)) /= size(var_perm,3) ) then
       call mpp_error(FATAL, "check_perm_3d: dimension mismatch")
    endif

    do k = 1, size(var,3)
       do j = 1, size(var,2)
          do i = 1, size(var,1)
             idx = [i,j,k]

             if (abs(var(i,j,k) - var_perm( idx(order(1)), idx(order(2)), idx(order(3)) ))  > 0) then

                print *, "perm mismatch at (x,y,z)=", i, j, k, "order=", order, &
                            " var  =", var(i,j,k), &
                            " perm =", var_perm(idx(order(1)), idx(order(2)), idx(order(3)))
                call mpp_error(FATAL, "check_perm_3d failed")
             endif
          enddo
       enddo
    enddo
  end subroutine check_perm_3d
end module
