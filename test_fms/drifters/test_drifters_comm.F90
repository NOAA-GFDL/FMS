
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

program test_drifters_comm
#ifdef use_drifters

  use drifters_core_mod
  use drifters_comm_mod
  use mpp_mod
  use mpp_domains_mod
  use fms_mod

  implicit none

  integer, parameter :: nd=2, npmax = 4
  integer :: nx, ny, halox, haloy, layout(2), i, j, npes, pe, it, nt
!  _TYPE_DOMAIN2D      :: domain
  type(domain2d)  :: domain
  type(drifters_core_type)      :: drfts
  type(drifters_comm_type) :: drfts_com
  real, parameter                     :: xmin=0., xmax=1., ymin=0., ymax=1.
  real                                :: dx, dy, u0, v0, dt, Lx, Ly
  real, allocatable                   :: x(:), y(:)
  character(len=128)  :: ermsg = ''
  integer :: ids(npmax)
  real    :: positions(nd, npmax), velocity(nd, npmax)
  integer :: io_status
!!$  integer :: stackmax=4000000

  namelist /drifters_comm_nml/ nx, ny, halox, haloy, u0, v0, dt, nt

  call mpp_init !(MPP_DEBUG)
  !call mpp_set_stack_size(3145746)

  ! default input values
  nx = 11
  ny = 21
  halox = 2
  haloy = 2
  u0    = 1.0
  v0    = 0.0
  dt    = 0.1
  nt    = 10

  ! read input
  read (input_nml_file, drifters_comm_nml, iostat=io_status)
  io_status = check_nml_error(io_status, 'test_drifters_comm')

  ! create global domain
  Lx = xmax - xmin
  Ly = ymax - ymin
  dx = Lx/real(nx-1)
  dy = Ly/real(ny-1)
  allocate(x(nx), y(ny))
  x = xmin + (/ ( i*dx, i=0, nx-1) /)
  y = ymin + (/ ( j*dy, j=0, ny-1) /)

  ! decompose domain
  call mpp_domains_init ! (MPP_DEBUG)
!!$  call mpp_domains_set_stack_size(stackmax)

  npes = mpp_npes()
  call mpp_define_layout( (/1,nx, 1,ny/), npes, layout )
  if(mpp_pe()==0) print *,'LAYOUT: ', layout
  call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halox, yhalo=haloy,&
       & xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN)

  ! set up drifters' communicator
  call drifters_comm_new(drfts_com)
  call drifters_comm_set_domain(drfts_com, domain, x, y, 0,0)
  ! repeated calls just for the fun of it
  call drifters_comm_set_domain(drfts_com, domain, x, y, 0,0)
  call drifters_comm_set_domain(drfts_com, domain, x, y, 0,0)

  ! create drifters
  call drifters_core_new(drfts, nd, npmax, ermsg)
  if(ermsg /= '') then
      print *,'ERROR',ermsg
      call mpp_error(FATAL, 'test_drifters_comm: Error')
  endif

  pe = mpp_pe()
  ids = (/ (i+100*pe, i=1,npmax) /)
  call drifters_core_set_ids(drfts, ids, ermsg)
  if(ermsg /= '') then
       print *,'ERROR', ermsg
       call mpp_error(FATAL, 'test_drifters_comm: Error')
  endif

  ! position particles
  if(pe == 0) then
     positions(:, 1) = (/ (drfts_com%xcmin + drfts_com%xcmax)/2., &
          &               (drfts_com%ycmin + drfts_com%ycmax)/2. /)
     !positions(:, 2) = (/0.,0.01/)
     call drifters_core_set_positions(drfts, positions(:, 1:1), ermsg)
     if(ermsg /= '') then
         print *,'ERROR',ermsg
         call mpp_error(FATAL, 'test_drifters_comm: Error')
     endif
  endif

  ! push drifters
  velocity(:,1) = (/u0, v0/)
  do it = 1, nt
     positions(:,1:drfts%np) = xmin + &
          & modulo(drfts%positions(:,1:drfts%np) + dt*velocity(:,1:drfts%np)-xmin, xmax-xmin)
     ! this will redistribute the drifters and update the positions
     call drifters_comm_update(drfts_com, drfts, positions(:,1:drfts%np))

     if(drfts%np > 0) then
        do i=1,drfts%np
           print '(a,i6,a,i3,a,i3,a, i3, a,2f10.6)', 'PE: ',pe, ' it=', it, ' np=', drfts%np, ' ip=', i, &
                & ' x,y=', drfts%positions(1,i), drfts%positions(2,i)
        enddo
     endif
!!$     call drifters_print(drfts, ermsg)
!!$     if(ermsg /= '') print *,ermsg
  enddo

  ! clean up
  call drifters_core_del(drfts, ermsg)
  if(ermsg /= '') then
       print *,ermsg, '<--- THIS'
       call mpp_error(FATAL, 'test_drifters_comm: Error')
  endif
  call drifters_comm_del(drfts_com)
  call mpp_domains_exit
  call mpp_exit

#endif
end program test_drifters_comm
