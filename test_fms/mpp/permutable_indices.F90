! This module allows arrays to be allocated and accessed using indices whose
! orders are allowed to differ from the index order of the underlying arrays.

! API:
! * allocate
! * section

module permutable_indices_mod
  implicit none

  type permutable_indices(ndim)
    integer, len :: ndim
    integer :: lb(ndim), ub(ndim)

    contains

    procedure :: permute => permutable_indices_permute
    procedure :: n => permutable_indices_n
  end type permutable_indices

  contains

  subroutine permutable_indices_permute(self, p)
    class(permutable_indices(*)), intent(inout) :: self
    integer, intent(in) :: p

    call permute_list(self%lb, p)
    call permute_list(self%ub, p)
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

  ! p should range from 1 to size(list)!
  subroutine permute_list(list, p)
    integer, intent(inout) :: list(:)
    integer, intent(in) :: p
    integer :: choices(size(list))
    integer :: n, k, i, f, indx

    n = size(list)
    if (p.lt.1 .or. p.gt.factorial(n)) then
      print *, "Error: p parameter is out of bounds"
      stop 1
    endif

    choices = list
    k = p - 1

    do i=1,n
      f = factorial(n - i)
      indx = k / f + 1
      k = mod(k, f)

      list(i) = choices(indx)
      choices(indx) = choices(n + 1 - i)
    enddo
  end subroutine permute_list
end module permutable_indices_mod
