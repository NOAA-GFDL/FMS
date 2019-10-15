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

program test_drifters

!* contents of input file: drifters_inp_test_3d.nc
!!$netcdf drifters_inp_test_3d {
!!$dimensions:
!!$   nd = 3 ; // number of dimensions (2 or 3)
!!$   np = 4 ; // number of particles
!!$variables:
!!$   double positions(np, nd) ;
!!$      positions:names = "x y z" ;
!!$      positions:units = "- - -" ;
!!$   int ids(np) ;
!!$
!!$// global attributes:
!!$      :velocity_names = "u v w" ;
!!$      :field_names = "temp" ;
!!$      :field_units = "C" ;
!!$      :time_units = "seconds" ;
!!$      :title = "example of input data for drifters" ;
!!$data:
!!$
!!$ positions =
!!$  -0.8, 0., 0.,
!!$  -0.2, 0., 0.,
!!$   0.2, 0., 0.,
!!$   0.8, 0., 0.;
!!$
!!$ ids = 1, 2, 3, 4 ; // must range from 1 to np, in any order
!!$}
!***********************************************************************

  ! Example showing how to use drifters_mod.

  use drifters_mod
  use mpp_mod
  use mpp_domains_mod

  implicit none

  ! declare drifters object
  type(drifters_type) :: drfts  ! drifters' object
  type(drifters_type) :: drfts2 ! to test copy
  character(len=128)  :: ermesg

  real    :: t0, dt, t, tend, rho
  real    :: xmin, xmax, ymin, ymax, zmin, zmax, theta
  real, parameter :: pi = 3.1415926535897931159980
  real, allocatable :: x(:), y(:)

!  if (DIMS == 2) then
!     real, allocatable :: u(:,:), v(:,:), w(:,:), temp(:,:)
!  endif

  real, allocatable :: z(:), u(:,:,:), v(:,:,:), w(:,:,:), temp(:,:,:)

  integer :: layout(2), nx, ny, nz, halox, haloy, i, j, k, npes, pe, root
  integer :: isd,  ied,  jsd,  jed, isc,  iec,  jsc,  jec
  integer :: pe_beg, pe_end
  integer :: ibnds(1) ! only used in _SERIAL mode

  type(domain2d) :: domain

  call mpp_init

  npes   = 1 !_MPP_NPES
  pe     = 0 !_MPP_PE
  root   = 0 !_MPP_ROOT
  pe_beg = npes/2
  pe_end = npes-1


  ! input parameters
  t0 = 0.0 ! initial time
  tend = 2.0*pi ! max time
  dt =  tend/20.0 ! time step
  ! domain boundaries
  xmin = -1. ; xmax = 1.
  ymin = -1. ; ymax = 1.
  zmin = -1. ; zmax = 1.
  nx = 41; ny = 41; nz = 21;
  halox = 2; haloy = 2;

  allocate( x(1-halox:nx+halox), y(1-haloy:ny+haloy))
  x = xmin + (xmax-xmin)*(/ (real(i-1)/real(nx-1), i = 1-halox, nx+halox) /)
  y = ymin + (ymax-ymin)*(/ (real(j-1)/real(ny-1), j = 1-haloy, ny+haloy) /)

!#if _DIMS == 2
!  allocate( u(1-halox:nx+halox, 1-haloy:ny+haloy), &
!       &    v(1-halox:nx+halox, 1-haloy:ny+haloy), &
!       &    w(1-halox:nx+halox, 1-haloy:ny+haloy), &
!       & temp(1-halox:nx+halox, 1-haloy:ny+haloy))
! #endif

!#if _DIMS == 3
  allocate( z(nz) )
  z = zmin + (zmax-zmin)*(/ (real(k-1)/real(nz-1), k = 1, nz) /)
  allocate( u(1-halox:nx+halox, 1-haloy:ny+haloy, nz), &
       &    v(1-halox:nx+halox, 1-haloy:ny+haloy, nz), &
       &    w(1-halox:nx+halox, 1-haloy:ny+haloy, nz), &
       & temp(1-halox:nx+halox, 1-haloy:ny+haloy, nz))
!#endif


!#ifndef _SERIAL
  ! decompose domain
  call mpp_domains_init ! (MPP_DEBUG)
!!$  call mpp_domains_set_stack_size(stackmax)

  call mpp_declare_pelist( (/ (i, i=pe_beg, pe_end) /), '_drifters')
!#endif

  ! this sumulates a run on a subset of PEs
  if(pe >= pe_beg .and. pe <= pe_end) then

!#ifndef _SERIAL
     call mpp_set_current_pelist( (/ (i, i=pe_beg, pe_end) /) )

     call mpp_define_layout( (/1,nx, 1,ny/), pe_end-pe_beg+1, layout )
     if(pe==root) print *,'LAYOUT: ', layout
     call mpp_define_domains((/1,nx, 1,ny/), layout, domain, &
          & xhalo=halox, yhalo=haloy) !,&
     !& xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN)
!#endif

     ! constructor
!#if _DIMS == 2
!     call drifters_new(drfts, &
!          & input_file ='drifters_inp_test_2d.nc'  , &
!          & output_file='drifters_out_test_2d.nc', &
!          & ermesg=ermesg)
!#endif

!#if _DIMS == 3
     call drifters_new(drfts, &
          & input_file ='drifters_inp_test_3d.nc'  , &
          & output_file='drifters_out_test_3d.nc', &
          & ermesg=ermesg)
!#endif
     if(ermesg/='') call my_error_handler(ermesg)

     ! set start/end pe
     drfts%comm%pe_beg = pe_beg
     drfts%comm%pe_end = pe_end

     ! set the initial time and dt
     drfts%time = t0
     drfts%dt   = dt

!#ifndef _SERIAL
     call mpp_get_data_domain   ( domain, isd,  ied,  jsd,  jed  )
     call mpp_get_compute_domain( domain, isc,  iec,  jsc,  jec  )
!#else
!     ibnds = lbound(x); isd = ibnds(1)
!     ibnds = ubound(x); ied = ibnds(1)
!     ibnds = lbound(y); jsd = ibnds(1)
!     ibnds = ubound(y); jed = ibnds(1)
!     isc = isd; iec = ied - 1
!     jsc = jsd; jec = jed - 1
!#endif


     ! set the PE domain boundaries. Xmin_comp/ymin_comp, xmax_comp/ymax_comp
     ! refer to the "compute" domain, which should cover densily the domain: ie
     ! xcmax[pe] = xcmin[pe_east]
     ! ycmax[pe] = ycmin[pe_north]
     ! Xmin_data/ymin_data, xmax_data/ymax_data refer to the "data" domain, which
     ! should be larger than the compute domain and therefore overlap: ie
     ! xdmax[pe] > xdmin[pe_east]
     ! ydmax[pe] > ydmin[pe_north]
     ! Particles in the overlap regions are tracked by several PEs.

     call drifters_set_domain(drfts, &
          & xmin_comp=x(isc  ), xmax_comp=x(iec+1), &
          & ymin_comp=y(jsc  ), ymax_comp=y(jec+1), &
          & xmin_data=x(isd  ), xmax_data=x(ied  ), &
          & ymin_data=y(jsd  ), ymax_data=y(jed  ), &
          !!$       & xmin_glob=xmin    , xmax_glob=xmax    , & ! periodicity in x
!!$       & ymin_glob=ymin    , ymax_glob=ymax    , & ! periodicity in y
     & ermesg=ermesg)
     if(ermesg/='') call my_error_handler(ermesg)

     ! set neighboring PEs [domain2d is of type(domain2d)]

     call drifters_set_pe_neighbors(drfts, domain=domain, ermesg=ermesg)
     if(ermesg/='') call my_error_handler(ermesg)

     ! set the velocities axes. Each velocity can have different axes.

     call drifters_set_v_axes(drfts, component='u', &
          & x=x, y=y, &
!#if _DIMS == 2
!          & z=DRFT_EMPTY_ARRAY, &
!#endif
!#if _DIMS >= 3
          & z=z, &
!#endif
          & ermesg=ermesg)
     if(ermesg/='') call my_error_handler(ermesg)

     call drifters_set_v_axes(drfts, component='v', &
          & x=x, y=y, &
!#if _DIMS == 2
!          & z=DRFT_EMPTY_ARRAY, &
!#endif

!#if _DIMS >= 3
          & z=z, &
!#endif
          & ermesg=ermesg)
     if(ermesg/='') call my_error_handler(ermesg)

!#if _DIMS == 3
     call drifters_set_v_axes(drfts, component='w', &
          & x=x, y=y, &
!#if _DIMS == 2
!          & z=DRFT_EMPTY_ARRAY, &
!#endif
!#if _DIMS >= 3
          & z=z, &
!#endif
          & ermesg=ermesg)
     if(ermesg/='') call my_error_handler(ermesg)
!#endif

     ! Distribute the drifters across PEs
     call drifters_distribute(drfts, ermesg)
     if(ermesg/='') call my_error_handler(ermesg)

     t = t0

     do while (t <= tend+epsilon(1.))

        ! Update time

        t = t + dt/2.0

        ! Set velocity and field
!#if _DIMS == 2
!        do j = 1-haloy, ny+haloy
!           do i = 1-halox, nx+halox
!              theta = atan2(y(j), x(i))
!              rho   = sqrt(x(i)**2 + y(j)**2)
!              u(i,j) = - rho * sin(theta)
!              v(i,j) = + rho * cos(theta)
!              temp(i,j) = (x(i)**2 + y(j)**2)
!           enddo
!        enddo
!        ! Push the drifters
!        call drifters_push(drfts, u=u, v=v, ermesg=ermesg)
!        if(ermesg/='') call my_error_handler(ermesg)

!#endif

!#if _DIMS == 3
        do k = 1, nz
           do j = 1-haloy, ny+haloy
              do i = 1-halox, nx+halox
                 theta = atan2(y(j), x(i))
                 rho   = sqrt(x(i)**2 + y(j)**2)
                 u(i,j,k) = - rho * sin(theta)
                 v(i,j,k) = + rho * cos(theta)
                 w(i,j,k) = + 0.01 * cos(t)
                 temp(i,j,k) = (x(i)**2 + y(j)**2) * (1.0 - z(k)**2)
              enddo
           enddo
        enddo
        ! Push the drifters
        call drifters_push(drfts, u=u, v=v, w=w, ermesg=ermesg)
        if(ermesg/='') call my_error_handler(ermesg)
!#endif


        ! Check if RK4 integration is complete

        if(drfts%rk4_completed) then

           ! Interpolate fields

           call drifters_set_field(drfts, index_field=1, x=x, y=y, &
!#if _DIMS >= 3
                & z=z, &
!#endif
                &    data=temp, ermesg=ermesg)
           if(ermesg/='') call my_error_handler(ermesg)

           ! Save data

           call drifters_save(drfts, ermesg=ermesg)
           if(ermesg/='') call my_error_handler(ermesg)

        endif

     enddo

     ! Write restart file

     call drifters_write_restart(drfts, filename='drifters_res.nc', &
          & ermesg=ermesg)
     if(ermesg/='') call my_error_handler(ermesg)

     ! test copy
     drfts2 = drfts

     ! destroy

     call drifters_del(drfts, ermesg=ermesg)
     if(ermesg/='') call my_error_handler(ermesg)

     deallocate(x, y)
     deallocate(u, v, temp)
!#if _DIMS == 3
     deallocate(z, w)
!#endif

  endif

!#ifndef _SERIAL
  call mpp_exit
!#endif

end program test_drifters

subroutine my_error_handler(mesg)
!#ifndef _SERIAL
  use mpp_mod, only : FATAL, mpp_error
!#endif
  implicit none
  character(len=*), intent(in) :: mesg
!#ifndef _SERIAL
  call mpp_error(FATAL, mesg)
!#else
!  print *, mesg
!  stop
!#endif

end subroutine my_error_handler
