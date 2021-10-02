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
!> @defgroup mpp_memutils_mod mpp_memutils_mod
!> @ingroup mpp
!> @brief Routines to initialize and report on memory usage during the model run.

!> @file
!> @brief File for mpp_memutils_mod

!> @addtogroup mpp_memutils_mod
!> @{
module mpp_memutils_mod

  use mpp_mod, only: mpp_min, mpp_max, mpp_sum, mpp_pe, mpp_root_pe
  use mpp_mod, only: mpp_error, FATAL, stderr, mpp_npes, get_unit
  use platform_mod

  implicit none
  private

  public :: mpp_print_memuse_stats, mpp_mem_dump
  public :: mpp_memuse_begin, mpp_memuse_end

  real    :: begin_memuse
  logical :: memuse_started = .false.

contains

  !#######################################################################
  !> Initialize the memory module, and record the initial memory use.
  subroutine mpp_memuse_begin
    if(memuse_started) then
       call mpp_error(FATAL, "mpp_memutils_mod: mpp_memuse_begin was already called")
    endif
    memuse_started = .true.

    call mpp_mem_dump(begin_memuse)
  end subroutine mpp_memuse_begin

  !#######################################################################
  !> End the memory collection, and report on total memory used during the
  !! execution of the model run.
  subroutine mpp_memuse_end(text, unit)
    character(len=*), intent(in) :: text !< Text to include in memory use statement
    integer, intent(in), optional :: unit !< Fortran unit number where memory report should go.
                                          !! Default is stderr.
    real    :: m, mmin, mmax, mavg, mstd, end_memuse
    integer :: mu

    if(.NOT.memuse_started) then
       call mpp_error(FATAL, "mpp_memutils_mod: mpp_memuse_begin must be called before calling mpp_memuse_being")
    endif
    memuse_started = .false.

    call mpp_mem_dump(end_memuse)

    mu = stderr(); if( PRESENT(unit) )mu = unit
    m = end_memuse - begin_memuse
    mmin = m; call mpp_min(mmin)
    mmax = m; call mpp_max(mmax)
    mavg = m; call mpp_sum(mavg); mavg = mavg/mpp_npes()
    mstd = (m-mavg)**2; call mpp_sum(mstd); mstd = sqrt( mstd/mpp_npes() )
    if( mpp_pe().EQ.mpp_root_pe() )write( mu,'(a64,4es11.3)' ) &
         'Memory(MB) used in '//trim(text)//'=', mmin, mmax, mstd, mavg

    return
  end subroutine mpp_memuse_end

  !#######################################################################
  !> Print the current memory high water mark to stderr, or the unit
  !! specified.
  subroutine mpp_print_memuse_stats(text, unit)
    character(len=*), intent(in) :: text !< Text to include in memory print statement
    integer, intent(in), optional :: unit !< Fortran unit number where print statement should go.
                                          !! Default is stderr.
    real :: m, mmin, mmax, mavg, mstd
    integer :: mu

    mu = stderr(); if( PRESENT(unit) )mu = unit
    call mpp_mem_dump(m)

    mmin = m; call mpp_min(mmin)
    mmax = m; call mpp_max(mmax)
    mavg = m; call mpp_sum(mavg); mavg = mavg/mpp_npes()
    mstd = (m-mavg)**2; call mpp_sum(mstd); mstd = sqrt( mstd/mpp_npes() )
    if( mpp_pe().EQ.mpp_root_pe() )write( mu,'(a64,4es11.3)' ) &
         'Memuse(MB) at '//trim(text)//'=', mmin, mmax, mstd, mavg

    return
  end subroutine mpp_print_memuse_stats

  !#######################################################################

  !> \brief Return the memory high water mark in MiB
  !!
  !! Query the system for the memory high water mark, return the result in MiB.
  subroutine mpp_mem_dump(memuse)
    real, intent(out) :: memuse !< Memory, high water mark, in MiB

    interface
      integer(KIND=c_size_t) function getpeakrss() bind(c, name="getpeakrss")
        use, intrinsic :: iso_c_binding
      end function getpeakrss
    end interface

    ! Value of Bytes to Mebibytes
    real, parameter :: B_to_MiB = 1048576.0

    ! Get the max memory use, convert to MiB
    memuse = real(getpeakrss())/B_to_MiB

    return
  end subroutine mpp_mem_dump
end module mpp_memutils_mod
!> @}
! close documentation grouping
