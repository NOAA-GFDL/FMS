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

!> @brief  This programs tests diag_manager with the following diag_table

program test_modern_diag
use   mpp_domains_mod,  only: domain2d, mpp_domains_set_stack_size, mpp_define_domains, mpp_define_io_domain, &
                              mpp_define_mosaic, domainug, mpp_get_compute_domains, mpp_define_unstruct_domain, &
                              mpp_get_compute_domain, mpp_get_data_domain, mpp_get_UG_domain_grid_index, &
                              mpp_get_UG_compute_domain
use   diag_manager_mod, only: diag_manager_init, diag_manager_end, diag_axis_init, register_diag_field, &
                              diag_axis_add_attribute, diag_field_add_attribute, diag_send_complete, &
                              diag_manager_set_time_end, send_data, register_static_field, &
                              diag_field_add_cell_measures
use   platform_mod,     only: r8_kind, r4_kind
use   fms_mod,          only: fms_init, fms_end
use   mpp_mod,          only: FATAL, mpp_error, mpp_npes, mpp_pe, mpp_root_pe, mpp_broadcast, input_nml_file
use   time_manager_mod, only: time_type, set_calendar_type, set_date, JULIAN, set_time, OPERATOR(+)
use fms_diag_object_mod,only: dump_diag_obj

implicit none

type(time_type)                   :: Time             !< Time of the simulation
type(time_type)                   :: Time_step        !< Time_step of the simulation
integer, dimension(2)             :: layout           !< Layout to use when setting up the domain
integer, dimension(2)             :: io_layout        !< io layout to use when setting up the io domain
integer                           :: nx               !< Number of x points
integer                           :: ny               !< Number of y points
integer                           :: nz               !< Number of z points
integer                           :: ug_dim_size      !< Number of points in the UG
type(domain2d)                    :: Domain           !< 2D domain
type(domain2d)                    :: Domain_cube_sph  !< cube sphere domain
type(domainug)                    :: land_domain      !< Unstructured domain
real,    dimension(:), allocatable:: x                !< X axis data
real,    dimension(:), allocatable:: y                !< Y axis_data
real,    dimension(:), allocatable:: z                !< Z axis_data
integer, dimension(:), allocatable:: ug_dim_data      !< UG axis_data
integer                           :: i                !< For do loops
integer                           :: id_x             !< axis id for the x dimension
integer                           :: id_y             !< axis id for the y dimension
integer                           :: id_UG            !< axis id for the unstructured dimension
integer                           :: id_z             !< axis id for the z dimention
integer                           :: id_lon
integer                           :: id_lat

integer                           :: id_var
real, dimension(:,:), allocatable :: var_data         !< Dummy variable data to send to diag_manager
logical                           :: used             !< Used for send_data call

call fms_init
call set_calendar_type(JULIAN)
call diag_manager_init

nx = 96
ny = 96
nz = 5
layout = (/1, mpp_npes()/)
io_layout = (/1, 1/)

!> Set up a normal (lat/lon) 2D domain, a cube sphere, and UG domain
call set_up_2D_domain(domain, layout, nx, ny, io_layout)
call set_up_cube_sph_domain(Domain_cube_sph, nx, ny, io_layout)
call create_land_domain(Domain_cube_sph, nx, ny, 6, land_domain, npes_group=1)
call mpp_get_UG_compute_domain(land_domain, size=ug_dim_size)

allocate(ug_dim_data(ug_dim_size))
call mpp_get_UG_domain_grid_index(land_domain, ug_dim_data)
ug_dim_data = ug_dim_data - 1

! Set up the data
allocate(x(ug_dim_size), y(ug_dim_size), z(nz))
do i=1,ug_dim_size
  x(i) = i
enddo
do i=1,ug_dim_size
  y(i) = i
enddo
do i=1,nz
   z(i) = i
enddo

! Set up the initial time
Time = set_date(2,1,1,0,0,0)

! Register the diags axis
id_x  = diag_axis_init('grid_xt',  x,  'point_E', 'x', long_name='point_E', set_name="land")
id_y  = diag_axis_init('grid_yt',  y,  'point_N', 'y', long_name='point_N', set_name="land")
id_z  = diag_axis_init('z',  z,  'point_Z', 'z', long_name='point_Z')

id_ug = diag_axis_init("grid_index",  real(ug_dim_data), "none", "U", long_name="grid indices", &
                         DomainU=land_domain, aux="geolon_t geolat_t", set_name="land")

call diag_axis_add_attribute (id_ug, 'compress', 'grid_xt grid_yt')

! Register the variables
id_lon = register_diag_field  ('lnd_mod', 'grid_xt', (/id_x/), Time, 'lon', 'mullions')
id_lat = register_diag_field  ('lnd_mod', 'grid_yt', (/id_y/), Time, 'lat', 'mullions')


id_var = register_diag_field  ('lnd_mod', 'var1', (/id_ug, id_z /), Time, 'Some scalar var', 'mullions', &
                                 standard_name="Land is important!")

call diag_manager_set_time_end(Time)
call diag_manager_set_time_end(set_date(2,1,2,0,0,0))

allocate(var_data(ug_dim_size, nz))

Time_step = set_time (3600,0) !< 1 hour

used = send_data(id_lon, x, Time)
used = send_data(id_lat, y, Time)

do i=1,23
  Time = Time + Time_step
  var_data = real(i)

  used = send_data(id_var, var_data, Time)

  call diag_send_complete(Time_step)
enddo

call diag_manager_end(Time)
call fms_end

contains

include "../fms2_io/create_atmosphere_domain.inc"
include "../fms2_io/create_land_domain.inc"

!> @brief Sets up a lat/lon domain
subroutine set_up_2D_domain(Domain, layout, nx, ny, io_layout)
  type(domain2d), intent(out)              :: Domain       !< 2D domain
  integer,        intent(in)               :: layout(:)    !< Layout to use when setting up the domain
  integer,        intent(in)               :: nx           !< Number of x points
  integer,        intent(in)               :: ny           !< Number of y points
  integer,        intent(in)               :: io_layout(:) !< Io layout to use when setting up the io_domain

  call mpp_domains_set_stack_size(17280000)
  call mpp_define_domains( (/1,nx,1,ny/), layout, Domain, name='2D domain')
  call mpp_define_io_domain(Domain, io_layout)
end subroutine set_up_2D_domain

!> @brief Sets up a cube sphere domain
subroutine set_up_cube_sph_domain(Domain_cube_sph, nx, ny, io_layout)
  type(domain2d), intent(out)              :: Domain_cube_sph  !< 2D domain
  integer,        intent(in)               :: nx               !< Number of x points
  integer,        intent(in)               :: ny               !< Number of y points
  integer,        intent(in)               :: io_layout(:)     !< Io layout to use when setting up the io_domain

  integer                                :: i              !< For do loops
  integer                                :: npes           !< Number of pes
  integer, parameter                     :: ntiles=6       !< Number of tiles
  integer,           dimension(4,ntiles) :: global_indices !< The global indices of each tile
  integer,           dimension(2,ntiles) :: layout         !< The layout of each tile
  integer,           dimension(ntiles)   :: pe_start       !< The starting PE of each tile
  integer,           dimension(ntiles)   :: pe_end         !< The ending PE of eeach tile

  npes = mpp_npes()

  !< Create the domain
  do i = 1,ntiles
    global_indices(:, i) = (/1, ny, 1, ny/)
    layout(:, i) = (/1, npes/ntiles/)
    pe_start(i) = (i-1)*(npes/ntiles)
    pe_end(i) = i*(npes/ntiles) - 1
  end do

  call create_atmosphere_domain((/nx, nx, nx, nx, nx, nx/), &
                                (/ny, ny, ny, ny, ny, ny/), &
                                global_indices, layout, pe_start, pe_end, &
                                io_layout, Domain_cube_sph)
end subroutine set_up_cube_sph_domain
end program test_modern_diag
