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

program test_modern_diag
#ifdef use_yaml
use   mpp_domains_mod,  only: domain2d, mpp_domains_set_stack_size, mpp_define_domains, mpp_define_io_domain, &
                              mpp_define_mosaic, domainug, mpp_get_compute_domains, mpp_define_unstruct_domain, &
                              mpp_get_compute_domain, mpp_get_data_domain, mpp_get_UG_domain_grid_index, &
                              mpp_get_UG_compute_domain
use   diag_manager_mod, only: diag_manager_init, diag_manager_end, diag_axis_init, register_diag_field, &
                              diag_axis_add_attribute, diag_field_add_attribute
use   fms_mod,          only: fms_init, fms_end
use   mpp_mod,          only: FATAL, mpp_error, mpp_npes, mpp_pe, mpp_root_pe, mpp_broadcast
use   time_manager_mod, only: time_type, set_calendar_type, set_date, JULIAN, set_time

implicit none

type(time_type)                   :: Time             !< Time of the simulation
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
integer                           :: id_x3            !< axis id for the x dimension in the cube sphere domain
integer                           :: id_y             !< axis id for the y dimension
integer                           :: id_y3            !< axis id for the y dimension in the cube sphere domain
integer                           :: id_UG            !< axis id for the unstructured dimension
integer                           :: id_z             !< axis id for the z dimention
integer                           :: id_var1          !< diag_field id for var in lon/lat grid
integer                           :: id_var2          !< diag_field id for var in lat/lon grid
integer                           :: id_var3          !< diag_field id for var in cube sphere grid
integer                           :: id_var4          !< diag_field id for 3d var in cube sphere grid
integer                           :: id_var5          !< diag_field id for var in UG grid
integer                           :: id_var6          !< diag_field id for var that is not domain decomposed
integer                           :: id_var7          !< Scalar var

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

! Set up the data
allocate(x(nx), y(ny), z(nz))
do i=1,nx
  x(i) = i
enddo
do i=1,ny
  y(i) = i
enddo
do i=1,nz
   z(i) = i
enddo

allocate(ug_dim_data(ug_dim_size))
call mpp_get_UG_domain_grid_index(land_domain, ug_dim_data)
ug_dim_data = ug_dim_data - 1

! Set up the intial time
Time = set_date(2,1,1,0,0,0)

! Register the diags axis
id_x  = diag_axis_init('x',  x,  'point_E', 'x', long_name='point_E', Domain2=Domain)
id_y  = diag_axis_init('y',  y,  'point_N', 'y', long_name='point_N', Domain2=Domain)

id_x3 = diag_axis_init('x3', x, 'point_E', 'x', Domain2=Domain_cube_sph)
id_y3 = diag_axis_init('y3', y, 'point_E', 'y', Domain2=Domain_cube_sph)

id_ug = diag_axis_init("grid_index",  real(ug_dim_data), "none", "U", long_name="grid indices", &
                         set_name="land", DomainU=land_domain, aux="geolon_t geolat_t")

id_z  = diag_axis_init('z',  z,  'point_Z', 'z', long_name='point_Z')
call diag_axis_add_attribute (id_z, 'formula', 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)')
call diag_axis_add_attribute (id_z, 'integer', 10)
call diag_axis_add_attribute (id_z, '1d integer', (/10, 10/))
call diag_axis_add_attribute (id_z, 'real', 10.)
call diag_axis_add_attribute (id_x, '1d real', (/10./))

if (id_x  .ne. 1) call mpp_error(FATAL, "The x axis does not have the expected id")
if (id_y  .ne. 2) call mpp_error(FATAL, "The y axis does not have the expected id")
if (id_x3 .ne. 3) call mpp_error(FATAL, "The x3 axis does not have the expected id")
if (id_y3 .ne. 4) call mpp_error(FATAL, "The y3 axis does not have the expected id")
if (id_ug .ne. 5) call mpp_error(FATAL, "The ug axis does not have the expected id")
if (id_z  .ne. 6) call mpp_error(FATAL, "The z axis does not have the expected id")

! Register the variables
id_var1 = register_diag_field  ('ocn_mod', 'var1', (/id_x, id_y/), Time, 'Var in a lon/lat domain', 'mullions')
id_var2 = register_diag_field  ('ocn_mod', 'var2', (/id_y, id_x/), Time,  &
                                'Var in a lon/lat domain with flipped dimensions', 'mullions')
id_var3 = register_diag_field  ('atm_mod', 'var3', (/id_x3, id_y3/), Time, 'Var in a cube sphere domain', 'mullions')
id_var4 = register_diag_field  ('atm_mod', 'var4', (/id_x3, id_y3, id_z/), Time, &
                                '3D var in a cube sphere domain', 'mullions')
id_var5 = register_diag_field  ('lnd_mod', 'var5', (/id_ug/), Time, 'Var in a UG domain', 'mullions')
id_var6 = register_diag_field  ('atm_mod', 'var6', (/id_z/), Time, 'Var not domain decomposed', 'mullions')

!< This has the same name as var1, but it should have a different id because the module is different
!! so it should have its own diag_obj
id_var7 = register_diag_field  ('lnd_mod', 'var1', Time, 'Some scalar var', 'mullions')

if (id_var1  .ne. 1) call mpp_error(FATAL, "var1 does not have the expected id")
if (id_var2  .ne. 2) call mpp_error(FATAL, "var2 does not have the expected id")
if (id_var3  .ne. 3) call mpp_error(FATAL, "var3 does not have the expected id")
if (id_var4  .ne. 4) call mpp_error(FATAL, "var4 does not have the expected id")
if (id_var5  .ne. 5) call mpp_error(FATAL, "var5 does not have the expected id")
if (id_var6  .ne. 6) call mpp_error(FATAL, "var6 does not have the expected id")
if (id_var7  .ne. 7) call mpp_error(FATAL, "var7 does not have the expected id")

call diag_field_add_attribute (id_var1, "some string", "this is a string")
call diag_field_add_attribute (id_var1, "integer", 10)
call diag_field_add_attribute (id_var1, "1d integer", (/10, 10/))
call diag_field_add_attribute (id_var1, "real", 10.)
call diag_field_add_attribute (id_var2, '1d real', (/10./))

call diag_manager_end(Time)
call fms_end

contains

include "../fms2_io/create_atmosphere_domain.inc"
include "../fms2_io/create_land_domain.inc"

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
#endif
end program test_modern_diag
