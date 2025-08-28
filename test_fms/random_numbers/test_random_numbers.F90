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
program test_random_numbers

use fms_mod, only: fms_init, fms_end
use platform_mod, only: r4_kind, r8_kind
use mpp_mod, only: mpp_error, fatal, stderr, mpp_npes, mpp_pe
use fms_string_utils_mod, only: string
use time_manager_mod, only: time_type, get_date, set_date, set_calendar_type, JULIAN
use random_numbers_mod

implicit none

integer, parameter :: n_moments = 1000 !> Highest-order moment to test
integer, parameter :: n0_1d = 2500 !> Initial length of 1D random sample vector
integer, parameter :: n0_2d = 50 !> Initial dimensions of 2D random sample array

integer, parameter :: seeds(*) = [0, -5, 3] !> Seed constants
integer, parameter :: k = TEST_FMS_KIND_ !> Either r4_kind or r8_kind

real(k), dimension(n_moments) :: moment_mu !> Expected moment values
real(k), dimension(n_moments) :: moment_sigma !> Standard deviations of sample moments

type(time_type) :: now !> Current date and time

call fms_init

if (mpp_npes() .ne. 4) then
  call mpp_error(fatal, "The random_numbers unit test requires four PEs")
endif

call set_time
call test_constructSeed
call test_getRandomNumbers

call fms_end

contains

!> Set `now` to the current date and time
subroutine set_time
  integer :: now_values(8) !> Values returned by the date_and_time() intrinsic

  call set_calendar_type(JULIAN)
  call date_and_time(values=now_values)

  now = set_date(now_values(1), now_values(2), now_values(3), &
                 now_values(5), now_values(6), now_values(7))
end subroutine set_time

!> Test random_numbers_mod's constructSeed()
subroutine test_constructSeed
  if (mpp_pe() .eq. 0) then
    call constructSeed_assert(seeds(1), seeds(1), now)
    call constructSeed_assert(seeds(1), seeds(1), now, 0)

    call constructSeed_assert(seeds(2), seeds(3), now)
    call constructSeed_assert(seeds(2), seeds(3), now, 0)

    call constructSeed_assert(seeds(1), seeds(1), now,  10)
    call constructSeed_assert(seeds(1), seeds(1), now, -10)

    call constructSeed_assert(seeds(2), seeds(3), now,  10)
    call constructSeed_assert(seeds(2), seeds(3), now, -10)

    call constructSeed_assert(seeds(1), seeds(1), now,  20)
    call constructSeed_assert(seeds(1), seeds(1), now, -20)

    call constructSeed_assert(seeds(2), seeds(3), now,  20)
    call constructSeed_assert(seeds(2), seeds(3), now, -20)
  endif
end subroutine test_constructSeed

subroutine constructSeed_assert(i, j, time, perm)
  integer,           intent(in) :: i, j !> Seed values
  type(time_type),   intent(in) :: time !> Time to be used in the construction of a seed vector
  integer, optional, intent(in) :: perm !> Permutation to be applied to the seed vector
  integer, dimension(8) :: seed_expected, seed_ret !> Expected and returned seed vectors
  integer :: year, month, day, hour, minute, second !> Date and time values

  call get_date(time, year, month, day, hour, minute, second)
  seed_expected = [i, j, year, month, day, hour, minute, second]
  if(present(perm)) seed_expected = ishftc(seed_expected, perm)

  seed_ret = constructSeed(i, j, time, perm)
  call array_compare_1d(seed_ret, seed_expected, "constructSeed unit test failed")
end subroutine constructSeed_assert

!> Test random_numbers_mod's getRandomNumbers()
subroutine test_getRandomNumbers
  type(randomNumberStream) :: stream !> Random number stream

  call calc_expected_moments

  select case (mpp_pe())
    case (0)
      stream = initializeRandomNumberStream(seeds(1))
      call test_getRandomNumbers_dispatch(stream)
    case (1)
      stream = initializeRandomNumberStream(seeds(2))
      call test_getRandomNumbers_dispatch(stream)
    case (2)
      stream = initializeRandomNumberStream(seeds)
      call test_getRandomNumbers_dispatch(stream)
    case (3)
      stream = initializeRandomNumberStream(constructSeed(seeds(2), seeds(3), now))
      call test_getRandomNumbers_dispatch(stream)
    case default
      call mpp_error(fatal, "Unexpected PE: " // string(mpp_pe()))
  end select
end subroutine test_getRandomNumbers

! Expression for the expectation value of the i-th raw moment
#define CALC_MOMENT_(i) (1._k / real(i + 1, k))

!> Calculate all expected moments and standard deviations of sample moments
subroutine calc_expected_moments
  integer :: i !> Moment order
  integer, parameter :: n = n0_1d !> Sample size for which to calculate moment standard deviations

  do i=1, n_moments
    moment_mu(i) = CALC_MOMENT_(i)
    moment_sigma(i) = sqrt((CALC_MOMENT_(2*i) - moment_mu(i)**2) / real(n, k))
  enddo
end subroutine calc_expected_moments

!> Invoke the 1D and 2D tests for a given random number stream
subroutine test_getRandomNumbers_dispatch(stream)
  type(randomNumberStream), intent(inout) :: stream !> Random number stream

  call test_samples_iter(stream, test_sample_1d, n0_1d)
  call test_samples_iter(stream, test_sample_2d, n0_2d)
end subroutine test_getRandomNumbers_dispatch

!> Run the requested test using progressively larger sample sizes, until ten
!> passing samples have been drawn
subroutine test_samples_iter(stream, test, n0)
  abstract interface
    !> Abstract interface for test_sample_1d and test_sample_2d
    function test_sample(stream, n)
      import :: randomNumberStream
      type(randomNumberStream), intent(inout) :: stream
      integer, intent(in) :: n
      logical :: test_sample
    end function
  end interface

  integer, parameter :: required_passes = 10 !> Number of samples that must pass the test

  type(randomNumberStream), intent(inout) :: stream !> Random number stream
  procedure(test_sample) :: test !> Function which draws a random sample and tests it
  integer, intent(in) :: n0 !> Initial sample size

  real(k) :: x !> Sample for 0D test
  integer :: n !> Sample size for 1D or 2D test
  integer :: pass_counter !> Number of test passes

  ! 0D case
  ! Draw a scalar and check that it's within [0,1]

  call getRandomNumbers(stream, x)
  call check_bounds(x)

  ! 1D and 2D cases
  ! Attempt to draw ten samples for which the first 1,000 moments are within one
  ! standard deviation of their expected values for the uniform distribution on
  ! [0,1]. Draw progressively larger samples until this condition has been fulfilled.

  n = n0
  pass_counter = 0

  do while (pass_counter .lt. required_passes)
    if (test(stream, n)) then
      pass_counter = pass_counter + 1
    endif

    n = n * 11 / 10
  enddo
end subroutine test_samples_iter

!> Draw a random sample and test its values and moments (1D)
function test_sample_1d(stream, n)
  type(randomNumberStream), intent(inout) :: stream !> Random number stream
  integer, intent(in) :: n !> Length of the 1D sample vector
  logical :: test_sample_1d !> True if the test passes for this sample
  real(k) :: v(n) !> Sample vector
  integer :: i !> Indices into v(:)

  call getRandomNumbers(stream, v)

  do i = 1, n
    call check_bounds(v(i))
  enddo

  test_sample_1d = compare_sample_moments(v)
end function test_sample_1d

!> Draw a random sample and test its values and moments (2D)
function test_sample_2d(stream, n)
  type(randomNumberStream), intent(inout) :: stream !> Random number stream
  integer, intent(in) :: n !> Dimensions of the 2D sample array
  logical :: test_sample_2d !> True if the test passes for this sample
  real(k) :: arr(n,n) !> Sample array
  integer :: i, j !> Indices into arr(:,:)

  call getRandomNumbers(stream, arr)

  do j = 1, n
    do i = 1, n
      call check_bounds(arr(i, j))
    enddo
  enddo

  test_sample_2d = compare_sample_moments(reshape(arr, [size(arr)]))
end function test_sample_2d

!> Check that the first 1,000 moments of a sample are within one standard
!> deviation of their expected values
function compare_sample_moments(v)
  real(k), intent(in) :: v(:) !> Vector containing a random sample
  logical :: compare_sample_moments !> True if the sample passes the test

  real(k), allocatable :: vi(:) !> v(:) raised to the power i
  integer :: i !> Moment order
  integer :: n !> Size of v(:)

  real(k) :: moment_sample !> Value of the i-th sample moment, calculated from v(:)

  n = size(v)
  allocate(vi(n))
  vi = 1._k

  do i = 1, n_moments
    vi = vi * v
    moment_sample = sum(vi) / n

    if (abs(moment_sample - moment_mu(i)) .gt. moment_sigma(i)) then
      compare_sample_moments = .false.
      return
    endif
  enddo

  compare_sample_moments = .true.
end function compare_sample_moments

!> Check that a value lies within the 0<x<1 range
subroutine check_bounds(x)
  real(k), intent(in) :: x !> Scalar value to check

  if (x.lt.0. .or. x.gt.1.) then
    call mpp_error(fatal, "Random number " // string(x) // " is out of expected bounds: [0, 1]")
  endif
end subroutine check_bounds

!> Compare two arrays of integers
subroutine array_compare_1d(arr1, arr2, msg)
  integer, intent(in), dimension(:) :: arr1, arr2 !> Arrays to be compared
  character(*), intent(in) :: msg !> Error message to be shown if the comparison fails
  integer :: m, n !> The sizes of the two arrays
  integer :: i !> Loop counter

  m = size(arr1)
  n = size(arr2)

  if (m .ne. n) then
    write(stderr(), "(A)") "1D array comparison failed due to incompatible array sizes"
    write(stderr(), "(A)") "Array 1 has size " // string(m) // " and array 2 has size " // string(n)
    call mpp_error(FATAL, msg)
  endif

  do i=1, m
    if (arr1(i) .ne. arr2(i)) then
      write(stderr(), "(A)") "1D array comparison failed due to element " // string(i)
      write(stderr(), "(A)") "Array 1 has value " // string(arr1(i)) // &
                           & " and array 2 has value " // string(arr2(i))
      call mpp_error(FATAL, msg)
    endif
  enddo
end subroutine array_compare_1d

end program test_random_numbers
