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

!> @brief Type to hold all the dummy data variables
type data_type
  real(kind=r8_kind), allocatable   :: var1(:,:)        !< Dummy data for var1
  real(kind=r8_kind), allocatable   :: var2(:,:)        !< Dummy data for var2
  real(kind=r8_kind), allocatable   :: var3(:,:)        !< Dummy data for var3
  real(kind=r8_kind), allocatable   :: var4(:,:,:)      !< Dummy data for var4
  real(kind=r8_kind), allocatable   :: var5(:)          !< Dummy data for var5
  real(kind=r8_kind), allocatable   :: var6(:)          !< Dummy data for var6
end type data_type

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
integer                           :: id_x3            !< axis id for the x dimension in the cube sphere domain
integer                           :: id_y             !< axis id for the y dimension
integer                           :: id_y3            !< axis id for the y dimension in the cube sphere domain
integer                           :: id_UG            !< axis id for the unstructured dimension
integer                           :: id_z             !< axis id for the z dimention
integer                           :: id_z2            !< axis id for the z dimention
integer                           :: id_var1          !< diag_field id for var in lon/lat grid
integer                           :: id_var2          !< diag_field id for var in lat/lon grid
integer                           :: id_var3          !< diag_field id for var in cube sphere grid
integer                           :: id_var4          !< diag_field id for 3d var in cube sphere grid
integer                           :: id_var5          !< diag_field id for var in UG grid
integer                           :: id_var6          !< diag_field id for var that is not domain decomposed
integer                           :: id_var7          !< 1D var
integer                           :: id_var8          !< Scalar var
type(data_type)                   :: var_data         !< Dummy variable data to send to diag_manager
logical                           :: used             !< Used for send_data call
integer                           :: io_status        !< Status after reading the namelist
logical                           :: debug = .false.  !< Flag used to ignore the axis/field_ids checks in the test.
                                                      !! Useful when using a portion or a different diag_table.yaml

namelist / test_modern_diag_nml / debug

call fms_init
call set_calendar_type(JULIAN)
call diag_manager_init

read (input_nml_file, test_modern_diag_nml, iostat=io_status)
if (io_status > 0) call mpp_error(FATAL,'=>test_modern_diag: Error reading input.nml')

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

id_z2 = diag_axis_init('z_edge',  z,  'point_Z', 'z', long_name='point_Z')
id_z  = diag_axis_init('z',  z,  'point_Z', 'z', long_name='point_Z', edges = id_z2)

call diag_axis_add_attribute (id_z, 'formula', 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)')
call diag_axis_add_attribute (id_z, 'integer', 10)
call diag_axis_add_attribute (id_z, '1d integer', (/10, 10/))
call diag_axis_add_attribute (id_z, 'real', 10.)
call diag_axis_add_attribute (id_x, '1d real', (/10./))
call diag_axis_add_attribute (id_ug, 'compress', 'x y')

if (.not. debug) then
  if (id_x  .ne. 1) call mpp_error(FATAL, "The x axis does not have the expected id")
  if (id_y  .ne. 2) call mpp_error(FATAL, "The y axis does not have the expected id")
  if (id_x3 .ne. 3) call mpp_error(FATAL, "The x3 axis does not have the expected id")
  if (id_y3 .ne. 4) call mpp_error(FATAL, "The y3 axis does not have the expected id")
  if (id_ug .ne. 5) call mpp_error(FATAL, "The ug axis does not have the expected id")
  if (id_z2 .ne. 6) call mpp_error(FATAL, "The z2 axis does not have the expected id")
  if (id_z  .ne. 7) call mpp_error(FATAL, "The z axis does not have the expected id")
endif

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
id_var8 = register_static_field ('atm_mod', 'var7', (/id_z/), "Be static!", "none")

if (.not. debug) then
  if (id_var1  .ne. 1) call mpp_error(FATAL, "var1 does not have the expected id")
  if (id_var2  .ne. 2) call mpp_error(FATAL, "var2 does not have the expected id")
  if (id_var3  .ne. 3) call mpp_error(FATAL, "var3 does not have the expected id")
  if (id_var4  .ne. 4) call mpp_error(FATAL, "var4 does not have the expected id")
  if (id_var5  .ne. 5) call mpp_error(FATAL, "var5 does not have the expected id")
  if (id_var6  .ne. 6) call mpp_error(FATAL, "var6 does not have the expected id")
  if (id_var7  .ne. 7) call mpp_error(FATAL, "var7 does not have the expected id")
  if (id_var8  .ne. 8) call mpp_error(FATAL, "var8 does not have the expected id")
endif

call diag_field_add_cell_measures(id_var6, area=id_var8, volume=id_var8)

call diag_field_add_attribute (id_var1, "some string", "this is a string")
call diag_field_add_attribute (id_var1, "integer", 10)
call diag_field_add_attribute (id_var1, "1d_integer", (/10, 10/))
call diag_field_add_attribute (id_var1, "real", 10.)
call diag_field_add_attribute (id_var2, '1d_real', (/10./))
call diag_field_add_attribute (id_var2, 'formula', 'p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)')
call diag_field_add_attribute (id_var2, 'cell_methods', 'area: mullions')

!! test dump routines
!! prints fields from objects for debugging to log if name is provided, othwerise goes to stdout
call dump_diag_obj('diag_obj_dump.log')
call dump_diag_obj()

call diag_manager_set_time_end(Time)
call diag_manager_set_time_end(set_date(2,1,2,0,0,0))

call allocate_dummy_data(var_data, domain, Domain_cube_sph, land_domain, nz)
Time_step = set_time (3600,0) !< 1 hour
do i=1,23
  Time = Time + Time_step
  call set_dummy_data(var_data, i)
  used = send_data(id_var1, var_data%var1, Time)
  used = send_data(id_var2, var_data%var2, Time)
  used = send_data(id_var3, var_data%var3, Time)
  used = send_data(id_var4, var_data%var4, Time)
  used = send_data(id_var5, var_data%var5, Time)
  used = send_data(id_var6, var_data%var6, Time)
  used = send_data(id_var7, var_data%var6, Time)

  !TODO I don't know about this (scalar field) or how this is suppose to work #WUT
  used = send_data(id_var8, var_data%var6, Time)

  call diag_send_complete(Time_step)
enddo
call deallocate_dummy_data(var_data)

call diag_manager_end(Time)
call fms_end

contains

include "../fms2_io/create_atmosphere_domain.inc"
include "../fms2_io/create_land_domain.inc"

!> @brief Allocates the dummy data to send to send_data
subroutine allocate_dummy_data(var, lat_lon_domain, cube_sphere, lnd_domain, nz)
  type(data_type), intent(inout) :: var             !< Data var to allocate
  type(domain2d),  intent(in)    :: lat_lon_domain  !< Lat/Lon domain
  type(domain2d),  intent(in)    :: cube_sphere     !< Cube sphere domain
  type(domainug),  intent(in)    :: lnd_domain      !< Land domain
  integer,         intent(in)    :: nz              !< Number of Z points

  integer :: nland !< Size of the unstructured grid per PE
  integer :: is    !< Starting x compute index
  integer :: ie    !< Ending x compute index
  integer :: js    !< Starting y compute index
  integer :: je    !< Ending y compute index

  call mpp_get_compute_domain(lat_lon_domain, is, ie, js, je)
  allocate(var%var1(is:ie, js:je)) !< Variable in a lat/lon domain
  allocate(var%var2(js:je, is:ie)) !< Variable in a lat/lon domain with flipped dimensions

  call mpp_get_compute_domain(cube_sphere, is, ie, js, je)
  allocate(var%var3(is:ie, js:je)) !< Variable in a cube sphere domain
  allocate(var%var4(is:ie, js:je, nz)) !< Variable in a 3D cube sphere domain

  call mpp_get_UG_compute_domain(lnd_domain, size=nland)
  allocate(var%var5(nland)) !< Variable in the land unstructured domain

  allocate(var%var6(nz)) !< 1D variable not domain decomposed

end subroutine allocate_dummy_data

!> @brief Allocates the dummy data to send to send_data
subroutine deallocate_dummy_data(var)
  type(data_type), intent(inout) :: var             !< Data var to deallocate

  deallocate(var%var1, var%var2, var%var3, var%var4, var%var5, var%var6)
end subroutine deallocate_dummy_data

!> @brief Sets the dummy_data to use in send_data
subroutine set_dummy_data(var, data_value)
  type(data_type), intent(inout) :: var        !< Data type to set
  integer,         intent(in)    :: data_value !< Value to send the data as

  var%var1 = real(data_value, kind=r8_kind)
  var%var2 = real(data_value + 1, kind=r8_kind)
  var%var3 = real(data_value + 2, kind=r8_kind)
  var%var4 = real(data_value + 3, kind=r8_kind)
  var%var5 = real(data_value + 4, kind=r8_kind)
  var%var6 = real(data_value + 5, kind=r8_kind)

end subroutine set_dummy_data

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
