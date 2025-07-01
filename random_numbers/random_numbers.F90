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

!> @defgroup random_numbers_mod random_numbers_mod
!> @ingroup random_numbers
!> @brief Generic module to wrap random number generators.
!!
!> The module defines a type that identifies the particular stream of random
!!  numbers, and has procedures for initializing it and getting real numbers
!!  in the range 0 to 1.
!!  This version uses the Mersenne Twister to generate random numbers on [0, 1].

module random_numbers_mod

  use MersenneTwister_mod, only: randomNumberSequence, & ! The random number engine.
                                 new_RandomNumberSequence, getRandomReal
  use time_manager_mod, only: time_type, get_date
  use platform_mod, only: r4_kind, r8_kind

  implicit none
  private

  !> @brief Type to hold a stream of randomly generated numbers
  !> @ingroup random_numbers_mod
  type randomNumberStream
    type(randomNumberSequence) :: theNumbers
  end type randomNumberStream

  !> Returns scalar, 1 or 2 D random real numbers
  !!
  !> @param stream @ref randomNumberStream to generate from
  !> @param[out] number output number(s)
  !> @ingroup random_numbers_mod
  interface getRandomNumbers
    module procedure :: get_random_number_0d_r4, get_random_number_0d_r8
    module procedure :: get_random_number_1d_r4, get_random_number_1d_r8
    module procedure :: get_random_number_2d_r4, get_random_number_2d_r8
  end interface getRandomNumbers

  !> Initializes stream for generating random numbers.
  !> @ingroup random_numbers_mod
  interface initializeRandomNumberStream
    module procedure initializeRandomNumberStream_S, initializeRandomNumberStream_V
  end interface initializeRandomNumberStream

  public :: randomNumberStream,                             &
            initializeRandomNumberStream, getRandomNumbers, &
            constructSeed

!> @addtogroup random_numbers_mod
!> @{

contains

  !> Initialization
  function initializeRandomNumberStream_S(seed) result(new)
    integer, intent( in)     :: seed
    type(randomNumberStream) :: new

    new%theNumbers = new_RandomNumberSequence(seed)

  end function initializeRandomNumberStream_S

  function initializeRandomNumberStream_V(seed) result(new)
    integer, dimension(:), intent( in) :: seed
    type(randomNumberStream)           :: new

    new%theNumbers = new_RandomNumberSequence(seed)
  end function initializeRandomNumberStream_V

  !> Constructs a unique seed from grid cell index and model date/time
  !!   The perm is supplied we generate a different seed by
  !!   circularly shifting the bits of the seed - this is useful
  !!   if we want to create more than one seed for a given
  !!   column and model date/time.
  !!   Note that abs(perm) must be <= the number of bits used
  !!   to represent the default integer (likely 32)
  function constructSeed(i, j, time, perm) result(seed)
    integer,           intent( in)  :: i, j
    type(time_type),   intent( in) :: time
    integer, optional, intent( in) :: perm
    integer, dimension(8) :: seed

    ! Local variables
    integer :: year, month, day, hour, minute, second


    call get_date(time, year, month, day, hour, minute, second)
    seed = (/ i, j, year, month, day, hour, minute, second /)
    if(present(perm)) seed = ishftc(seed, perm)
  end function constructSeed

#include "random_numbers_r4.fh"
#include "random_numbers_r8.fh"

end module random_numbers_mod

!> @}
! close documentation grouping
