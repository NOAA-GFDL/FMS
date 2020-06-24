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

!> This is a simple test of MPP layouts and a simple domain. In
!> testing it should be run with 4 pes.

!> @author Ed Hartnett 6/22/20
program test_domains_simple
  use mpp_mod
  use mpp_domains_mod

  implicit none
#include "../../include/fms_platform.h"
  integer :: pe, npes
  integer :: nx=40, ny=40, nz=40
  integer :: layout(2)
  type(domain2D) :: domain_2D
  type(domain1D) :: domain_1D
  integer :: is, ie, js, je
  
  call mpp_init()

  pe = mpp_pe()
  npes = mpp_npes()

  !--- initialize mpp domains
  call mpp_domains_init(MPP_DEBUG)

  ! Define a layout.
  call mpp_define_layout((/1, nx, 1, ny/), npes, layout)

  ! Check results when run on 4 pes.
  if (npes .eq. 4) then
     if (layout(1) .ne. 2 .or. layout(2) .ne. 2) call mpp_error(FATAL, 'bad layout values')
  endif

  ! Define a 1D domain.
  !call mpp_define_domains((/1, nx/), 4, domain_1D)
  call mpp_define_domains((/1, nx/), 4, domain_1D, pelist=(/0, 1, 2, 3/))
  call mpp_get_compute_domain(domain_1D, is, ie)
  
  ! Check results when run on 4 pes.
  if (npes .eq. 4) then
     if (is .ne. pe * 10 + 1) call mpp_error(FATAL, 'bad domain start value')
     if (ie .ne. 10 * (pe + 1)) call mpp_error(FATAL, 'bad domain end value')
  endif

  ! Define a 2D domain.
  call mpp_define_domains((/1, nx, 1, ny/), layout, domain_2D)
  call mpp_get_compute_domain(domain_2D, xbegin = is, xend = ie, ybegin = js, yend = je)
  print *, pe, is, ie, js, je

  ! Check results when run on 4 pes.
  if (npes .eq. 4) then
     if (pe .eq. 0 .or. pe .eq. 2) then
        if (is .ne. 1 .or. ie .ne. 20) call mpp_error(FATAL, 'bad 2D domain')
     else
       if (is .ne. 21 .or. ie .ne. 40) call mpp_error(FATAL, 'bad 2D domain')
     endif
     if (pe .eq. 0 .or. pe .eq. 1) then
        if (js .ne. 1 .or. je .ne. 20) call mpp_error(FATAL, 'bad 2D domain')
     else
        if (js .ne. 21 .or. je .ne. 40) call mpp_error(FATAL, 'bad 2D domain')
     endif
  endif

  ! Clean up domains.
  call mpp_domains_exit()

  ! Finalize MPP.
  call mpp_exit()

end program test_domains_simple
