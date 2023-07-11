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

!> @brief  This programs tests the modern diag_manager

program test_diag_openmp
  use omp_lib
  use mpp_mod,           only: mpp_npes, mpp_pe, mpp_sync
  use platform_mod,      only: r8_kind
  use mpp_domains_mod,   only: domain2d, mpp_define_domains, mpp_define_io_domain, mpp_get_compute_domain
  use block_control_mod, only: block_control_type, define_blocks
  use fms_mod,           only: fms_init, fms_end
  use diag_manager_mod,  only: diag_manager_init, diag_manager_end, diag_axis_init, register_diag_field, &
                               diag_send_complete, diag_manager_set_time_end, send_data, register_static_field
  use time_manager_mod,  only: time_type, set_calendar_type, set_date, JULIAN, set_time


  implicit none

  integer                            :: nx           !< Number of points in the x direction
  integer                            :: ny           !< Number of points in the y direction
  integer                            :: nz           !< Number of points in the z direction
  integer                            :: layout(2)    !< Layout
  integer                            :: io_layout(2) !< Io layout
  type(domain2d)                     :: Domain       !< 2D domain
  integer                            :: is           !< Starting x compute index
  integer                            :: ie           !< Ending x compute index
  integer                            :: js           !< Starting y compute index
  integer                            :: je           !< Ending y compute index
  type(time_type)                    :: Time         !< Time of the simulation
  type(time_type)                    :: Time_step    !< Time of the simulation
  real,    dimension(:), allocatable :: x            !< X axis data
  integer                            :: id_x         !< axis id for the x dimension
  real,    dimension(:), allocatable :: y            !< Y axis_data
  integer                            :: id_y         !< axis id for the y dimension
  real,    dimension(:), allocatable :: z            !< Z axis data
  integer                            :: id_z         !< axis id for the z dimension
  real(kind=r8_kind),    allocatable :: var(:,:,:)   !< Dummy variable data
  integer                            :: i, j         !< For do loops
  type(block_control_type)           :: my_block     !< Returns instantiated @ref block_control_type
  logical                            :: message      !< Flag for outputting debug message
  integer                            :: isw          !< Starting index for each thread in the x direction
  integer                            :: iew          !< Ending index for each thread in the x direction
  integer                            :: jsw          !< Starting index for each thread in the y direction
  integer                            :: jew          !< Ending index for each thread in the y direction
  integer                            :: id_var1      !< diag_field id for var in 1d
  integer                            :: id_var2      !< diag_field id for var in lon/lat grid
  integer                            :: id_var3      !< diag_field id for var in lon/lat/z grid
  logical                            :: used         !< .true. if the send_data call was sucessful

  call fms_init
  call set_calendar_type(JULIAN)
  call diag_manager_init

  nx = 96
  ny = 96
  nz = 5
  layout = (/1, mpp_npes()/)
  io_layout = (/1, 1/)

  ! Set up the intial time
  Time = set_date(2,1,1,0,0,0)

  !< Create a lat/lon domain
  call mpp_define_domains( (/1,nx,1,ny/), layout, Domain, name='2D domain')
  call mpp_define_io_domain(Domain, io_layout)
  call mpp_get_compute_domain(Domain, is, ie, js, je)

  ! Set up the data
  allocate(x(nx), y(ny), z(nz))
  allocate(var(is:ie, js:je, nz))
  do i=1,nx
    x(i) = i
  enddo

  do i=1,ny
    y(i) = i
  enddo

  do i=1,nz
    z(i) = i
  enddo

  !< Register the axis:
  id_x  = diag_axis_init('x',  x,  'point_E', 'x', long_name='point_E', Domain2=Domain)
  id_y  = diag_axis_init('y',  y,  'point_N', 'y', long_name='point_N', Domain2=Domain)
  id_z  = diag_axis_init('z',  z,  'pressure', 'z', long_name='too much pressure')

  !< Register the variables
  id_var1 = register_diag_field  ('ocn_mod', 'var1', (/id_x/), Time, 'Var in a lon domain', 'mullions')
  id_var2 = register_diag_field  ('ocn_mod', 'var2', (/id_x, id_y/), Time, 'Var in a lon/lat domain', 'mullions')
  id_var3 = register_diag_field  ('ocn_mod', 'var3', (/id_x, id_y, id_z/), Time, &
    'Var in a lon/lat/z domain', 'mullions')

  call diag_manager_set_time_end(set_date(2,1,2,0,0,0))

  !< Divide the domain further into blocks
  call define_blocks ('testing_model', my_block, is, ie, js, je, kpts=0, &
                      nx_block=1, ny_block=4, message=message)

  Time_step = set_time (3600,0) !< 1 hour
  do j = 1, 23 !simulated time
    Time = set_date(2,1,1,j,0,0)
    var = real(j, kind=r8_kind) !< Set the data
!$OMP parallel do default(shared) private(i, isw, iew, jsw, jew) schedule (dynamic,1)
    do i = 1, 4
      isw = my_block%ibs(i)
      jsw = my_block%jbs(i)
      iew = my_block%ibe(i)
      jew = my_block%jbe(i)

      used=send_data(id_var1, var(isw:iew, 1, 1), time, is_in=isw, ie_in=iew)
      used=send_data(id_var2, var(isw:iew, jsw:jew, 1), time, is_in=isw, js_in=jsw, &
                     ie_in=iew, je_in=jew)
      used=send_data(id_var3, var(isw:iew, jsw:jew, :), time, is_in=isw, js_in=jsw, &
                     ie_in=iew, je_in=jew, ks_in=1, ke_in=nz)
    enddo
    call diag_send_complete(Time_step)
  enddo

  call diag_manager_end(Time)
  call fms_end
end program test_diag_openmp