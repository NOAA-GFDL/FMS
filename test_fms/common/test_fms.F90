! This module allows arrays to be permuted, and provides a data type for the
! purpose of storing permuted array bounds. It provides procedures for
! initializing a 2D or 3D array with random data, and for comparing a 2D or
! 3D array with reference answers.

module fms_test_mod
  use random_numbers_mod, only: randomNumberStream, initializeRandomNumberStream, getRandomNumbers
  use mpp_mod,         only: mpp_error, FATAL
  use platform_mod

  implicit none

  interface arr_init
    module procedure :: arr_init_2d_r4, arr_init_2d_r8, arr_init_2d_i4, arr_init_2d_i8
    module procedure :: arr_init_3d_r4, arr_init_3d_r8, arr_init_3d_i4, arr_init_3d_i8
  end interface arr_init

  interface arr_compare
    module procedure :: arr_compare_2d_r4, arr_compare_2d_r8, arr_compare_2d_i4, arr_compare_2d_i8
    module procedure :: arr_compare_3d_r4, arr_compare_3d_r8, arr_compare_3d_i4, arr_compare_3d_i8
    module procedure :: arr_compare_4d_r4, arr_compare_4d_r8, arr_compare_4d_i4, arr_compare_4d_i8
  end interface arr_compare

  interface arr_compare_tol
    module procedure :: arr_compare_tol_2d_r4, arr_compare_tol_2d_r8
    module procedure :: arr_compare_tol_3d_r4, arr_compare_tol_3d_r8
    module procedure :: arr_compare_tol_4d_r4, arr_compare_tol_4d_r8

    module procedure :: arr_compare_tol_2d_scalar_r4, arr_compare_tol_2d_scalar_r8
    module procedure :: arr_compare_tol_3d_scalar_r4, arr_compare_tol_3d_scalar_r8
    module procedure :: arr_compare_tol_4d_scalar_r4, arr_compare_tol_4d_scalar_r8
  end interface arr_compare_tol

  type permutable_indices(ndim)
    integer, len :: ndim
    integer :: lb(ndim), ub(ndim)

    contains

    procedure :: permute => permutable_indices_permute
    procedure :: n => permutable_indices_n
  end type permutable_indices

  contains

#define FMS_TEST_TYPE_ real
#define TYPECAST_ real

#define FMS_TEST_KIND_ r4_kind

#define ARR_INIT_2D_ arr_init_2d_r4
#define ARR_INIT_3D_ arr_init_3d_r4
#define ARR_COMPARE_2D_ arr_compare_2d_r4
#define ARR_COMPARE_3D_ arr_compare_3d_r4
#define ARR_COMPARE_4D_ arr_compare_4d_r4
#include "include/test_fms.inc"
#undef ARR_INIT_2D_
#undef ARR_INIT_3D_
#undef ARR_COMPARE_2D_
#undef ARR_COMPARE_3D_
#undef ARR_COMPARE_4D_

#define ARR_COMPARE_TOL_2D_ arr_compare_tol_2d_r4
#define ARR_COMPARE_TOL_3D_ arr_compare_tol_3d_r4
#define ARR_COMPARE_TOL_4D_ arr_compare_tol_4d_r4
#define ARR_COMPARE_TOL_2D_SCALAR_ arr_compare_tol_2d_scalar_r4
#define ARR_COMPARE_TOL_3D_SCALAR_ arr_compare_tol_3d_scalar_r4
#define ARR_COMPARE_TOL_4D_SCALAR_ arr_compare_tol_4d_scalar_r4
#include "include/test_fms_real.inc"
#undef ARR_COMPARE_TOL_2D_
#undef ARR_COMPARE_TOL_3D_
#undef ARR_COMPARE_TOL_4D_
#undef ARR_COMPARE_TOL_2D_SCALAR_
#undef ARR_COMPARE_TOL_3D_SCALAR_
#undef ARR_COMPARE_TOL_4D_SCALAR_

#undef FMS_TEST_KIND_
#define FMS_TEST_KIND_ r8_kind

#define ARR_INIT_2D_ arr_init_2d_r8
#define ARR_INIT_3D_ arr_init_3d_r8
#define ARR_COMPARE_2D_ arr_compare_2d_r8
#define ARR_COMPARE_3D_ arr_compare_3d_r8
#define ARR_COMPARE_4D_ arr_compare_4d_r8
#include "include/test_fms.inc"
#undef ARR_INIT_2D_
#undef ARR_INIT_3D_
#undef ARR_COMPARE_2D_
#undef ARR_COMPARE_3D_
#undef ARR_COMPARE_4D_
#undef ARR_COMPARE_TOL_2D_SCALAR_
#undef ARR_COMPARE_TOL_3D_SCALAR_
#undef ARR_COMPARE_TOL_4D_SCALAR_

#define ARR_COMPARE_TOL_2D_ arr_compare_tol_2d_r8
#define ARR_COMPARE_TOL_3D_ arr_compare_tol_3d_r8
#define ARR_COMPARE_TOL_4D_ arr_compare_tol_4d_r8
#define ARR_COMPARE_TOL_2D_SCALAR_ arr_compare_tol_2d_scalar_r8
#define ARR_COMPARE_TOL_3D_SCALAR_ arr_compare_tol_3d_scalar_r8
#define ARR_COMPARE_TOL_4D_SCALAR_ arr_compare_tol_4d_scalar_r8
#include "include/test_fms_real.inc"
#undef ARR_COMPARE_TOL_2D_
#undef ARR_COMPARE_TOL_3D_
#undef ARR_COMPARE_TOL_4D_
#undef ARR_COMPARE_TOL_2D_SCALAR_
#undef ARR_COMPARE_TOL_3D_SCALAR_
#undef ARR_COMPARE_TOL_4D_SCALAR_

#undef FMS_TEST_KIND_

#undef FMS_TEST_TYPE_
#undef TYPECAST_

#define FMS_TEST_TYPE_ integer
#define TYPECAST_ int

#define FMS_TEST_KIND_ i4_kind
#define ARR_INIT_2D_ arr_init_2d_i4
#define ARR_INIT_3D_ arr_init_3d_i4
#define ARR_COMPARE_2D_ arr_compare_2d_i4
#define ARR_COMPARE_3D_ arr_compare_3d_i4
#define ARR_COMPARE_4D_ arr_compare_4d_i4
#include "include/test_fms.inc"
#undef FMS_TEST_KIND_
#undef ARR_INIT_2D_
#undef ARR_INIT_3D_
#undef ARR_COMPARE_2D_
#undef ARR_COMPARE_3D_
#undef ARR_COMPARE_4D_

#define FMS_TEST_KIND_ i8_kind
#define ARR_INIT_2D_ arr_init_2d_i8
#define ARR_INIT_3D_ arr_init_3d_i8
#define ARR_COMPARE_2D_ arr_compare_2d_i8
#define ARR_COMPARE_3D_ arr_compare_3d_i8
#define ARR_COMPARE_4D_ arr_compare_4d_i8
#include "include/test_fms.inc"
#undef FMS_TEST_KIND_
#undef ARR_INIT_2D_
#undef ARR_INIT_3D_
#undef ARR_COMPARE_2D_
#undef ARR_COMPARE_3D_
#undef ARR_COMPARE_4D_

#undef FMS_TEST_TYPE_
#undef TYPECAST_

  subroutine permutable_indices_permute(self, p)
    class(permutable_indices(*)), intent(inout) :: self
    integer, intent(in) :: p

    call permute_arr(self%lb, p)
    call permute_arr(self%ub, p)
  end subroutine permutable_indices_permute

  function permutable_indices_n(self, i) result(n)
    class(permutable_indices(*)), intent(inout) :: self
    integer, intent(in) :: i
    integer :: n

    n = self%ub(i) - self%lb(i) + 1
  end function permutable_indices_n

  pure recursive function factorial(n) result(res)
    integer, intent(in) :: n
    integer :: res

    if (n.eq.0) then
      res = 1
    else
      res = n * factorial(n-1)
    endif
  end function factorial

  subroutine permute_arr(arr, p)
    integer, intent(inout) :: arr(:) !< List to be permuted
    integer, intent(in) :: p !< Which permutation to produce: may range from 1 to size(arr)!
    integer :: choices(size(arr))
    integer :: n, k, i, f, indx

    n = size(arr)
    if (p.lt.1 .or. p.gt.factorial(n)) then
      print *, "Error: p parameter is out of bounds"
      stop 1
    endif

    choices = arr
    k = p - 1

    do i=1,n
      f = factorial(n - i)
      indx = k / f + 1
      k = mod(k, f)

      arr(i) = choices(indx)
      choices(indx) = choices(n + 1 - i)
    enddo
  end subroutine permute_arr
end module fms_test_mod
