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

!> @file
!> This is a simple test of MPP layouts and a simple domain.

!> @author Ed Hartnett 6/22/20
program test_domains_simple
  use mpp_mod
  use mpp_domains_mod

  implicit none
#include "../../include/fms_platform.h"
  integer :: pe, npes
  integer :: nx=128, ny=128, nz=40, stackmax=4000000
  integer :: layout(2)
  type(domain2D) :: domain
  
  call mpp_init()

  pe = mpp_pe()
  npes = mpp_npes()

  !--- initialize mpp domains
  call mpp_domains_init(MPP_DEBUG)
  call mpp_domains_set_stack_size(stackmax)

  call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
  if (npes .eq. 4) then
     if (layout(1) .ne. 2 .or. layout(2) .ne. 2) call mpp_error(FATAL, 'bad layout values')
  endif
  
  call mpp_define_domains( (/1,nx,1,ny/), layout, domain )
  
  call mpp_domains_exit()
  call mpp_exit()

end program test_domains_simple
