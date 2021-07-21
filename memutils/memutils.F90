!***********************************************************************
!!                   GNU Lesser General Public License
!!
!! This file is part of the GFDL Flexible Modeling System (FMS).
!!
!! FMS is free software: you can redistribute it and/or modify it under
!! the terms of the GNU Lesser General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or (at
!! your option) any later version.
!!
!! FMS is distributed in the hope that it will be useful, but WITHOUT
!! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!! FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!! for more details.
!!
!! You should have received a copy of the GNU Lesser General Public
!! License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @defgroup memutils_mod memutils_mod
!> @ingroup memutils
!! @brief Module to expose the memory printing API
!! @author V. Balaji
!!
!! Module to expose the memory printing API
!!
!! Model components at times desire to print the memory statics to stderr,
!! or another file (e.g. the log file).  This can be useful in debugging.
!! This module exposes the print_memuse_stat and memutils_init calls for
!! use in external user code.

!> @file
!> @brief File for @ref memutils_mod

!> @addtogroup memutils_mod
!> @{
module memutils_mod
!Author: Balaji (V.Balaji@noaa.gov)
  use mpp_mod, only: mpp_pe, mpp_root_pe, mpp_npes, mpp_min, mpp_max, mpp_sum, stderr
  use mpp_memutils_mod, only: mpp_print_memuse_stats

  implicit none

  logical :: memutils_initialized=.FALSE. !< Indicate if module has been initialized.
  logical, private :: print_memory_usage=.FALSE. !< Default behavior of print_memuse_stats()

  public :: memutils_init
  public :: print_memuse_stats

contains

  !> @brief Initialize the memutils module
  !!
  !! memutils_init initializes the print_memory_usage and memutils_initialized
  !! module variables when called.  This configures if print_memuse_stats() should
  !! print the memory statistics when called.  memutils_init, at present, can
  !! be called multiple times, and potentially change the behavior of print_memuse_stats().
  subroutine memutils_init(print_flag)
    logical, optional :: print_flag !< Indicate if memory statistics should be printed by default

    if ( PRESENT(print_flag) ) print_memory_usage = print_flag
    memutils_initialized = .TRUE.
    return
  end subroutine memutils_init

!> @brief Print memory usage stats to stdout, or a particular file
!!
!! API to allow external user code to print memory statistics to stderr, or an
!! optional file (Fortran file unit).  The module variable `print_memory_usage`
!! will control the default behavior (set when memutils_init is called).  The
!! default behavior can be modified passing in the optional `always` (logical)
!! parameter.  If `always.eqv..TRUE.`, print_memuse_stats() will print the memory
!! statistics.
  subroutine print_memuse_stats( text, unit, always )
    character(len=*), intent(in) :: text !< Text to be printed before the memory statistics
    integer, intent(in), optional :: unit !< Fortran file unit to where memory statistics should be recorded
    logical, intent(in), optional :: always !< If `.TRUE.`, force memory statistics to be printed

    if( PRESENT(always) )then
      if( .NOT.always )return
    else
      if( .NOT.print_memory_usage )return
    end if
    call mpp_print_memuse_stats(text, unit)
    return
  end subroutine print_memuse_stats
end module memutils_mod
!> @}
! close documentation grouping
