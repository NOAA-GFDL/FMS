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
!> @defgroup mpp_efp_mod mpp_efp_mod
!> @ingroup mpp
!> @brief This module provides interfaces to the non-domain-oriented communication
!! subroutines.
!!
!> Mainly includes interfaces and type definitions for reproducing operations with extended
!! fixed point data.

!> @file
!> @brief file for @ref mpp_efp_mod

!> @addtogroup mpp_efp_mod
!> @{
module mpp_efp_mod

use mpp_mod, only : mpp_error, FATAL, WARNING, NOTE
use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_npes
use mpp_mod, only : mpp_sum
use platform_mod

implicit none ; private

public :: mpp_reproducing_sum, mpp_efp_list_sum_across_PEs
public :: mpp_efp_plus, mpp_efp_minus, mpp_efp_to_real, mpp_real_to_efp, mpp_efp_real_diff
public :: operator(+), operator(-), assignment(=)
public :: mpp_query_efp_overflow_error, mpp_reset_efp_overflow_error

integer, parameter :: NUMBIT = 46  !< number of bits used in the 64-bit signed integer representation.
integer, parameter :: NUMINT = 6   !< The number of long integers to use to represent
                                   !! a real number.

integer(i8_kind), parameter :: prec=2_8**NUMBIT !< The precision of each integer.
real(r8_kind), parameter :: r_prec=2.0_8**NUMBIT !< A real version of prec.
real(r8_kind), parameter :: I_prec=1.0_8/(2.0_8**NUMBIT) !< The inverse of prec.
integer, parameter :: max_count_prec=2**(63-NUMBIT)-1  !< The number of values that can be added together
                              !! with the current value of prec before there will
                              !! be roundoff problems.

real(r8_kind), parameter, dimension(NUMINT) :: &
  pr = (/ r_prec**2, r_prec, 1.0_8, 1.0_8/r_prec, 1.0_8/r_prec**2, 1.0_8/r_prec**3 /)
real(r8_kind), parameter, dimension(NUMINT) :: &
  I_pr = (/ 1.0_8/r_prec**2, 1.0_8/r_prec, 1.0_8, r_prec, r_prec**2, r_prec**3 /)

logical :: overflow_error = .false., NaN_error = .false.
logical :: debug = .false.    !< Making this true enables debugging output.

!> @}

!> This interface uses a conversion to an integer representation
!! of real numbers to give order-invariant sums that will reproduce
!! across PE count.
!!
!! This idea comes from R. Hallberg and A. Adcroft.
!!
!> @ingroup mpp_efp_mod
interface mpp_reproducing_sum
  module procedure mpp_reproducing_sum_r8_2d
  module procedure mpp_reproducing_sum_r8_3d
  module procedure mpp_reproducing_sum_r4_2d
end interface mpp_reproducing_sum

!> The Extended Fixed Point (mpp_efp) type provides a public interface for doing
!! sums and taking differences with this type.
!> @ingroup mpp_efp_mod
type, public :: mpp_efp_type
  private
  integer(i8_kind), dimension(NUMINT) :: v
end type mpp_efp_type


!> Operator override interface for mpp_efp_type
!> @ingroup mpp_efp_mod
interface operator (+); module procedure mpp_efp_plus  ; end interface
!> Operator override interface for mpp_efp_type
!> @ingroup mpp_efp_mod
interface operator (-); module procedure mpp_efp_minus ; end interface
!> Assignment override interface for mpp_efp_type
!> @ingroup mpp_efp_mod
interface assignment(=); module procedure mpp_efp_assign ; end interface

!> @addtogroup mpp_efp_mod
!> @{

contains

!> @brief Calculates a reproducing sum for a 2D, 8-byte real array
function mpp_reproducing_sum_r8_2d(array, isr, ier, jsr, jer, EFP_sum, reproducing, &
                            overflow_check, err) result(sum)
  real(r8_kind), dimension(:,:), intent(in) :: array
  integer,        optional,          intent(in) :: isr, ier, jsr, jer
  type(mpp_efp_type), optional,     intent(out) :: EFP_sum
  logical,        optional,          intent(in) :: reproducing
  logical,        optional,          intent(in) :: overflow_check
  integer,        optional,         intent(out) :: err
  real(r8_kind)                             :: sum  !< Result

  integer(i8_kind), dimension(NUMINT)  :: ints_sum
  integer(i8_kind) :: ival, prec_error
  real(r8_kind)  :: rsum(1), rs
  real(r8_kind)  :: max_mag_term
  logical :: repro, over_check
  character(len=256) :: mesg
  integer :: i, j, n, is, ie, js, je, sgn

  if (mpp_npes() > max_count_prec) call mpp_error(FATAL, &
    "mpp_reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/mpp_npes.")

  prec_error = (2_8**62 + (2_8**62 - 1)) / mpp_npes()

  is = 1 ; ie = size(array,1) ; js = 1 ; je = size(array,2 )
  if (present(isr)) then
    if (isr < is) call mpp_error(FATAL, &
      "Value of isr too small in mpp_reproducing_sum_2d.")
    is = isr
  endif
  if (present(ier)) then
    if (ier > ie) call mpp_error(FATAL, &
      "Value of ier too large in mpp_reproducing_sum_2d.")
    ie = ier
  endif
  if (present(jsr)) then
    if (jsr < js) call mpp_error(FATAL, &
      "Value of jsr too small in mpp_reproducing_sum_2d.")
    js = jsr
  endif
  if (present(jer)) then
    if (jer > je) call mpp_error(FATAL, &
      "Value of jer too large in mpp_reproducing_sum_2d.")
    je = jer
  endif

  repro = .true. ; if (present(reproducing)) repro = reproducing
  over_check = .true. ; if (present(overflow_check)) over_check = overflow_check

  if (repro) then
    overflow_error = .false. ; NaN_error = .false. ; max_mag_term = 0.0
    ints_sum(:) = 0
    if (over_check) then
      if ((je+1-js)*(ie+1-is) < max_count_prec) then
        do j=js,je ; do i=is,ie
          call increment_ints_faster(ints_sum, array(i,j), max_mag_term);
        enddo ; enddo
        call carry_overflow(ints_sum, prec_error)
      elseif ((ie+1-is) < max_count_prec) then
        do j=js,je
          do i=is,ie
            call increment_ints_faster(ints_sum, array(i,j), max_mag_term);
          enddo
          call carry_overflow(ints_sum, prec_error)
        enddo
      else
        do j=js,je ; do i=is,ie
          call increment_ints(ints_sum, real_to_ints(array(i,j), prec_error), &
                              prec_error);
        enddo ; enddo
      endif
    else
      do j=js,je ; do i=is,ie
        sgn = 1 ; if (array(i,j)<0.0) sgn = -1
        rs = abs(array(i,j))
        do n=1,NUMINT
          ival = int(rs*I_pr(n), 8)
          rs = rs - ival*pr(n)
          ints_sum(n) = ints_sum(n) + sgn*ival
        enddo
      enddo ; enddo
      call carry_overflow(ints_sum, prec_error)
    endif

    if (present(err)) then
      err = 0
      if (overflow_error) &
        err = err+2
      if (NaN_error) &
        err = err+4
      if (err > 0) then ; do n=1,NUMINT ; ints_sum(n) = 0 ; enddo ; endif
    else
      if (NaN_error) then
        call mpp_error(FATAL, "NaN in input field of mpp_reproducing_sum(_2d), this indicates numerical instability")
      endif
      if (abs(max_mag_term) >= prec_error*pr(1)) then
        write(mesg, '(ES13.5)') max_mag_term
        call mpp_error(FATAL,"Overflow in mpp_reproducing_sum(_2d) conversion of "//trim(mesg))
      endif
      if (overflow_error) then
        call mpp_error(FATAL, "Overflow in mpp_reproducing_sum(_2d).")
      endif
    endif

    call mpp_sum(ints_sum, NUMINT)

    call regularize_ints(ints_sum)
    sum = ints_to_real(ints_sum)
  else
    rsum(1) = 0.0
    do j=js,je ; do i=is,ie
      rsum(1) = rsum(1) + array(i,j);
    enddo ; enddo
    call mpp_sum(rsum,1)
    sum = rsum(1)

    if (present(err)) then ; err = 0 ; endif

    if (debug .or. present(EFP_sum)) then
      overflow_error = .false.
      ints_sum = real_to_ints(sum, prec_error, overflow_error)
      if (overflow_error) then
        if (present(err)) then
          err = err + 2
        else
          write(mesg, '(ES13.5)') sum
          call mpp_error(FATAL,"Repro_sum_2d: Overflow in real_to_ints conversion of "//trim(mesg))
        endif
      endif
    endif
  endif

  if (present(EFP_sum)) EFP_sum%v(:) = ints_sum(:)

  if (debug) then
    write(mesg,'("2d RS: ", ES24.16, 6 Z17.16)') sum, ints_sum(1:NUMINT)
    if(mpp_pe() == mpp_root_pe()) call mpp_error(NOTE, mesg)
  endif

end function mpp_reproducing_sum_r8_2d

function mpp_reproducing_sum_r4_2d(array, isr, ier, jsr, jer, EFP_sum, reproducing, &
                            overflow_check, err) result(sum)
  real(r4_kind), dimension(:,:), intent(in) :: array
  integer,        optional,          intent(in) :: isr, ier, jsr, jer
  type(mpp_efp_type), optional,     intent(out) :: EFP_sum
  logical,        optional,          intent(in) :: reproducing
  logical,        optional,          intent(in) :: overflow_check
  integer,        optional,         intent(out) :: err
  real(r4_kind)                             :: sum  !< Result

  real(r8_kind) :: array_r8(size(array,1), size(array,2))

  array_r8 = array

  sum = mpp_reproducing_sum_r8_2d(array_r8, isr, ier, jsr, jer, EFP_sum, reproducing, &
                            overflow_check, err)

  return

end function mpp_reproducing_sum_r4_2d

!> @brief Reproducing sum for 3d arrays of 8-bit reals
!!
!> This function uses a conversion to an integer representation
!! of real numbers to give order-invariant sums that will reproduce
!! across PE count.  This idea comes from R. Hallberg and A. Adcroft.
function mpp_reproducing_sum_r8_3d(array, isr, ier, jsr, jer, sums, EFP_sum, err) &
                            result(sum)
  real(r8_kind), dimension(:,:,:),        intent(in) :: array
  integer,    optional,                       intent(in) :: isr, ier, jsr, jer
  real(r8_kind), dimension(:), optional, intent(out) :: sums
  type(mpp_efp_type),     optional,          intent(out) :: EFP_sum
  integer,            optional,              intent(out) :: err
  real(r8_kind)                                      :: sum  !< Result

  real(r8_kind)    :: max_mag_term
  integer(i8_kind), dimension(NUMINT)  :: ints_sum
  integer(i8_kind), dimension(NUMINT,size(array,3))  :: ints_sums
  integer(i8_kind) :: prec_error
  character(len=256) :: mesg
  integer :: i, j, k, is, ie, js, je, ke, isz, jsz, n

  if (mpp_npes() > max_count_prec) call mpp_error(FATAL, &
    "mpp_reproducing_sum: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/mpp_npes.")

  prec_error = (2_8**62 + (2_8**62 - 1)) / mpp_npes()
  max_mag_term = 0.0

  is = 1 ; ie = size(array,1) ; js = 1 ; je = size(array,2) ; ke = size(array,3)
  if (present(isr)) then
    if (isr < is) call mpp_error(FATAL, &
      "Value of isr too small in mpp_reproducing_sum(_3d).")
    is = isr
  endif
  if (present(ier)) then
    if (ier > ie) call mpp_error(FATAL, &
      "Value of ier too large in mpp_reproducing_sum(_3d).")
    ie = ier
  endif
  if (present(jsr)) then
    if (jsr < js) call mpp_error(FATAL, &
      "Value of jsr too small in mpp_reproducing_sum(_3d).")
    js = jsr
  endif
  if (present(jer)) then
    if (jer > je) call mpp_error(FATAL, &
      "Value of jer too large in mpp_reproducing_sum(_3d).")
    je = jer
  endif
  jsz = je+1-js; isz = ie+1-is

  if (present(sums)) then
    if (size(sums) > ke) call mpp_error(FATAL, "Sums is smaller than "//&
      "the vertical extent of array in mpp_reproducing_sum(_3d).")
    ints_sums(:,:) = 0
    overflow_error = .false. ; NaN_error = .false. ; max_mag_term = 0.0
    if (jsz*isz < max_count_prec) then
      do k=1,ke
        do j=js,je ; do i=is,ie
          call increment_ints_faster(ints_sums(:,k), array(i,j,k), max_mag_term);
        enddo ; enddo
        call carry_overflow(ints_sums(:,k), prec_error)
      enddo
    elseif (isz < max_count_prec) then
      do k=1,ke ; do j=js,je
        do i=is,ie
          call increment_ints_faster(ints_sums(:,k), array(i,j,k), max_mag_term);
        enddo
        call carry_overflow(ints_sums(:,k), prec_error)
      enddo ; enddo
    else
      do k=1,ke ; do j=js,je ; do i=is,ie
        call increment_ints(ints_sums(:,k), &
                            real_to_ints(array(i,j,k), prec_error), prec_error);
      enddo ; enddo ; enddo
    endif
    if (present(err)) then
      err = 0
      if (abs(max_mag_term) >= prec_error*pr(1)) err = err+1
      if (overflow_error) err = err+2
      if (NaN_error) err = err+2
      if (err > 0) then ; do k=1,ke ; do n=1,NUMINT ; ints_sums(n,k) = 0 ; enddo ; enddo ; endif
    else
      if (NaN_error) call mpp_error(FATAL, &
             "NaN in input field of mpp_reproducing_sum(_3d), this indicates numerical instability")
      if (abs(max_mag_term) >= prec_error*pr(1)) then
        write(mesg, '(ES13.5)') max_mag_term
        call mpp_error(FATAL,"Overflow in mpp_reproducing_sum(_3d) conversion of "//trim(mesg))
      endif
      if (overflow_error) call mpp_error(FATAL, "Overflow in mpp_reproducing_sum(_3d).")
    endif

    call mpp_sum(ints_sums(:,1:ke), NUMINT*ke)

    sum = 0.0
    do k=1,ke
      call regularize_ints(ints_sums(:,k))
      sums(k) = ints_to_real(ints_sums(:,k))
      sum = sum + sums(k)
    enddo

    if (present(EFP_sum)) then
      EFP_sum%v(:) = 0
      do k=1,ke ; call increment_ints(EFP_sum%v(:), ints_sums(:,k)) ; enddo
    endif

    if (debug) then
      do n=1,NUMINT ; ints_sum(n) = 0 ; enddo
      do k=1,ke ; do n=1,NUMINT ; ints_sum(n) = ints_sum(n) + ints_sums(n,k) ; enddo ; enddo
      write(mesg,'("3D RS: ", ES24.16, 6 Z17.16)') sum, ints_sum(1:NUMINT)
      if(mpp_pe()==mpp_root_pe()) call mpp_error(NOTE, mesg)
    endif
  else
    ints_sum(:) = 0
    overflow_error = .false. ; NaN_error = .false. ; max_mag_term = 0.0
    if (jsz*isz < max_count_prec) then
      do k=1,ke
        do j=js,je ; do i=is,ie
          call increment_ints_faster(ints_sum, array(i,j,k), max_mag_term);
        enddo ; enddo
        call carry_overflow(ints_sum, prec_error)
      enddo
    elseif (isz < max_count_prec) then
      do k=1,ke ; do j=js,je
        do i=is,ie
          call increment_ints_faster(ints_sum, array(i,j,k), max_mag_term);
        enddo
        call carry_overflow(ints_sum, prec_error)
      enddo ; enddo
    else
      do k=1,ke ; do j=js,je ; do i=is,ie
        call increment_ints(ints_sum, real_to_ints(array(i,j,k), prec_error), &
                            prec_error);
      enddo ; enddo ; enddo
    endif
    if (present(err)) then
      err = 0
      if (abs(max_mag_term) >= prec_error*pr(1)) err = err+1
      if (overflow_error) err = err+2
      if (NaN_error) err = err+2
      if (err > 0) then ; do n=1,NUMINT ; ints_sum(n) = 0 ; enddo ; endif
    else
      if (NaN_error) call mpp_error(FATAL, &
          "NaN in input field of mpp_reproducing_sum(_3d), this indicates numerical instability")
      if (abs(max_mag_term) >= prec_error*pr(1)) then
        write(mesg, '(ES13.5)') max_mag_term
        call mpp_error(FATAL,"Overflow in mpp_reproducing_sum(_3d) conversion of "//trim(mesg))
      endif
      if (overflow_error) call mpp_error(FATAL, "Overflow in mpp_reproducing_sum(_3d).")
    endif

    call mpp_sum(ints_sum, NUMINT)

    call regularize_ints(ints_sum)
    sum = ints_to_real(ints_sum)

    if (present(EFP_sum)) EFP_sum%v(:) = ints_sum(:)

    if (debug) then
      write(mesg,'("3d RS: ", ES24.16, 6 Z17.16)') sum, ints_sum(1:NUMINT)
      if(mpp_pe()==mpp_root_pe()) call mpp_error(NOTE, mesg)
    endif
  endif

end function mpp_reproducing_sum_r8_3d

!> @brief This function converts a real number to an equivalent representation
!! using several long integers.
function real_to_ints(r, prec_error, overflow) result(ints)
  real(r8_kind),            intent(in) :: r
  integer(i8_kind), optional, intent(in) :: prec_error
  logical,   optional,       intent(inout) :: overflow
  integer(i8_kind),    dimension(NUMINT) :: ints

  real(r8_kind) :: rs
  character(len=80) :: mesg
  integer(i8_kind) :: ival, prec_err
  integer :: sgn, i

  prec_err = prec ; if (present(prec_error)) prec_err = prec_error
  ints(:) = 0_8
  if ((r >= 1e30) .eqv. (r < 1e30)) then ; NaN_error = .true. ; return ; endif

  sgn = 1 ; if (r<0.0) sgn = -1
  rs = abs(r)

  if (present(overflow)) then
    if (.not.(rs < prec_err*pr(1))) overflow = .true.
    if ((r >= 1e30) .eqv. (r < 1e30)) overflow = .true.
  elseif (.not.(rs < prec_err*pr(1))) then
    write(mesg, '(ES13.5)') r
    call mpp_error(FATAL,"Overflow in real_to_ints conversion of "//trim(mesg))
  endif

  do i=1,NUMINT
    ival = int(rs*I_pr(i), 8)
    rs = rs - ival*pr(i)
    ints(i) = sgn*ival
  enddo

end function real_to_ints

!> @brief This function reverses the conversion in real_to_ints.
function ints_to_real(ints) result(r)
  integer(i8_kind), dimension(NUMINT), intent(in) :: ints
  real(r8_kind) :: r

  integer :: i

  r = 0.0
  do i=1,NUMINT ; r = r + pr(i)*ints(i) ; enddo
end function ints_to_real

!> @brief This subroutine increments a number with another, both using the integer
!! representation in real_to_ints.
subroutine increment_ints(int_sum, int2, prec_error)
  integer(i8_kind), dimension(NUMINT), intent(inout) :: int_sum
  integer(i8_kind), dimension(NUMINT), intent(in)    :: int2
  integer(i8_kind), optional,      intent(in)    :: prec_error

  integer :: i

  do i=NUMINT,2,-1
    int_sum(i) = int_sum(i) + int2(i)
    ! Carry the local overflow.
    if (int_sum(i) > prec) then
      int_sum(i) = int_sum(i) - prec
      int_sum(i-1) = int_sum(i-1) + 1
    elseif (int_sum(i) < -prec) then
      int_sum(i) = int_sum(i) + prec
      int_sum(i-1) = int_sum(i-1) - 1
    endif
  enddo
  int_sum(1) = int_sum(1) + int2(1)
  if (present(prec_error)) then
    if (abs(int_sum(1)) > prec_error) overflow_error = .true.
  else
    if (abs(int_sum(1)) > prec) overflow_error = .true.
  endif

end subroutine increment_ints

!> @brief This subroutine increments a number with another, both using the integer
!! representation in real_to_ints, but without doing any carrying of overflow.
!! The entire operation is embedded in a single call for greater speed.
subroutine increment_ints_faster(int_sum, r, max_mag_term)
  integer(i8_kind), dimension(NUMINT), intent(inout) :: int_sum
  real(r8_kind),                        intent(in) :: r
  real(r8_kind),                     intent(inout) :: max_mag_term

  real(r8_kind) :: rs
  integer(i8_kind) :: ival
  integer :: sgn, i

  if ((r >= 1e30) .eqv. (r < 1e30)) then ; NaN_error = .true. ; return ; endif
  sgn = 1 ; if (r<0.0) sgn = -1
  rs = abs(r)
  if (rs > abs(max_mag_term)) max_mag_term = r

  do i=1,NUMINT
    ival = int(rs*I_pr(i), 8)
    rs = rs - ival*pr(i)
    int_sum(i) = int_sum(i) + sgn*ival
  enddo

end subroutine increment_ints_faster

!> @brief This subroutine handles carrying of the overflow.
subroutine carry_overflow(int_sum, prec_error)
  integer(i8_kind), dimension(NUMINT), intent(inout) :: int_sum
  integer(i8_kind),                intent(in)    :: prec_error

  integer :: i, num_carry

  do i=NUMINT,2,-1 ; if (abs(int_sum(i)) > prec) then
    num_carry = int(int_sum(i) * I_prec)
    int_sum(i) = int_sum(i) - num_carry*prec
    int_sum(i-1) = int_sum(i-1) + num_carry
  endif ; enddo
  if (abs(int_sum(1)) > prec_error) then
    overflow_error = .true.
  endif

end subroutine carry_overflow

!> @brief This subroutine carries the overflow, and then makes sure that
!! all integers are of the same sign as the overall value.
subroutine regularize_ints(int_sum)
  integer(i8_kind), dimension(NUMINT), intent(inout) :: int_sum

  logical :: positive
  integer :: i, num_carry

  do i=NUMINT,2,-1 ; if (abs(int_sum(i)) > prec) then
    num_carry = int(int_sum(i) * I_prec)
    int_sum(i) = int_sum(i) - num_carry*prec
    int_sum(i-1) = int_sum(i-1) + num_carry
  endif ; enddo

  ! Determine the sign of the final number.
  positive = .true.
  do i=1,NUMINT
    if (abs(int_sum(i)) > 0) then
      if (int_sum(i) < 0) positive = .false.
      exit
    endif
  enddo

  if (positive) then
    do i=NUMINT,2,-1 ; if (int_sum(i) < 0) then
      int_sum(i) = int_sum(i) + prec
      int_sum(i-1) = int_sum(i-1) - 1
    endif ; enddo
  else
    do i=NUMINT,2,-1 ; if (int_sum(i) > 0) then
      int_sum(i) = int_sum(i) - prec
      int_sum(i-1) = int_sum(i-1) + 1
    endif ; enddo
  endif

end subroutine regularize_ints

function mpp_query_efp_overflow_error()
  logical :: mpp_query_efp_overflow_error
  mpp_query_efp_overflow_error = overflow_error
end function mpp_query_efp_overflow_error

subroutine mpp_reset_efp_overflow_error()
  overflow_error = .false.
end subroutine mpp_reset_efp_overflow_error

function mpp_efp_plus(EFP1, EFP2)
  type(mpp_efp_type)             :: mpp_efp_plus
  type(mpp_efp_type), intent(in) :: EFP1, EFP2

  mpp_efp_plus = EFP1

  call increment_ints(mpp_efp_plus%v(:), EFP2%v(:))
end function mpp_efp_plus

function mpp_efp_minus(EFP1, EFP2)
  type(mpp_efp_type)             :: mpp_efp_minus
  type(mpp_efp_type), intent(in) :: EFP1, EFP2
  integer :: i

  do i=1,NUMINT ; mpp_efp_minus%v(i) = -1*EFP2%v(i) ; enddo

  call increment_ints(mpp_efp_minus%v(:), EFP1%v(:))
end function mpp_efp_minus

!> @brief This subroutine assigns all components of the extended fixed point type
!! variable on the RHS (EFP2) to the components of the variable on the LHS
!! (EFP1).
subroutine mpp_efp_assign(EFP1, EFP2)
  type(mpp_efp_type), intent(out) :: EFP1
  type(mpp_efp_type), intent(in)  :: EFP2
  integer i

  do i=1,NUMINT ; EFP1%v(i) = EFP2%v(i) ; enddo
end subroutine mpp_efp_assign

function mpp_efp_to_real(EFP1)
  type(mpp_efp_type), intent(inout) :: EFP1
  real(r8_kind) :: mpp_efp_to_real

  call regularize_ints(EFP1%v)
  mpp_efp_to_real = ints_to_real(EFP1%v)
end function mpp_efp_to_real

function mpp_efp_real_diff(EFP1, EFP2)
  type(mpp_efp_type), intent(in) :: EFP1, EFP2
  real(r8_kind) :: mpp_efp_real_diff

  type(mpp_efp_type)             :: EFP_diff

  EFP_diff = EFP1 - EFP2
  mpp_efp_real_diff = mpp_efp_to_real(EFP_diff)

end function mpp_efp_real_diff

function mpp_real_to_efp(val, overflow)
  real(r8_kind),    intent(in) :: val
  logical, optional, intent(inout) :: overflow
  type(mpp_efp_type)               :: mpp_real_to_efp

  logical :: over
  character(len=80) :: mesg

  if (present(overflow)) then
    mpp_real_to_efp%v(:) = real_to_ints(val, overflow=overflow)
  else
    over = .false.
    mpp_real_to_efp%v(:) = real_to_ints(val, overflow=over)
    if (over) then
      write(mesg, '(ES13.5)') val
      call mpp_error(FATAL,"Overflow in mpp_real_to_efp conversion of "//trim(mesg))
    endif
  endif

end function mpp_real_to_efp

!> This subroutine does a sum across PEs of a list of EFP variables,
!! returning the sums in place, with all overflows carried.
subroutine mpp_efp_list_sum_across_PEs(EFPs, nval, errors)
  type(mpp_efp_type), dimension(:), intent(inout) :: EFPs
  integer, intent(in) :: nval
  logical, dimension(:), optional, intent(out) :: errors

  integer(i8_kind), dimension(NUMINT,nval) :: ints
  integer(i8_kind) :: prec_error
  logical :: error_found
  character(len=256) :: mesg
  integer :: i, n

  if (mpp_npes() > max_count_prec) call mpp_error(FATAL, &
    "mpp_efp_list_sum_across_PEs: Too many processors are being used for the value of "//&
    "prec.  Reduce prec to (2^63-1)/mpp_npes.")

  prec_error = (2_8**62 + (2_8**62 - 1)) / mpp_npes()
  ! overflow_error is an overflow error flag for the whole module.
  overflow_error = .false. ; error_found = .false.

  do i=1,nval ; do n=1,NUMINT ; ints(n,i) = EFPs(i)%v(n) ; enddo ; enddo

  call mpp_sum(ints(:,:), NUMINT*nval)

  if (present(errors)) errors(:) = .false.
  do i=1,nval
    overflow_error = .false.
    call carry_overflow(ints(:,i), prec_error)
    do n=1,NUMINT ; EFPs(i)%v(n) = ints(n,i) ; enddo
    if (present(errors)) errors(i) = overflow_error
    if (overflow_error) then
      write (mesg,'("mpp_efp_list_sum_across_PEs error at ",i6," val was ",ES12.6, ", prec_error = ",ES12.6)') &
             i, mpp_efp_to_real(EFPs(i)), real(prec_error)
      if(mpp_pe()==mpp_root_pe()) call mpp_error(WARNING, mesg)
    endif
    error_found = error_found .or. overflow_error
  enddo
  if (error_found .and. .not.(present(errors))) then
    call mpp_error(FATAL, "Overflow in mpp_efp_list_sum_across_PEs.")
  endif

end subroutine mpp_efp_list_sum_across_PEs

end module mpp_efp_mod
!> @}
! close documentation grouping
