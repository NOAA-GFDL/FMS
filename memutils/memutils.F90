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
module memutils_mod
!Author: Balaji (V.Balaji@noaa.gov)
  use mpp_mod, only: mpp_pe, mpp_root_pe, mpp_npes, mpp_min, mpp_max, mpp_sum, stderr
  use mpp_memutils_mod, only: mpp_print_memuse_stats

  logical :: memutils_initialized=.FALSE.

  public :: memutils_init
  public :: print_memuse_stats

  logical, private :: print_memory_usage=.FALSE.

contains

  !> \brief Initialize the memutils module
  subroutine memutils_init(print_flag)
    logical, optional :: print_flag

    if ( PRESENT(print_flag) ) print_memory_usage = print_flag
    memutils_initialized = .TRUE.
    return
  end subroutine memutils_init

  subroutine print_memuse_stats( text, unit, always )
    character(len=*), intent(in) :: text
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: always

    if( PRESENT(always) )then
      if( .NOT.always )return
    else
      if( .NOT.print_memory_usage )return
    end if
    call mpp_print_memuse_stats(text, unit)
    return
  end subroutine print_memuse_stats
end module memutils_mod
