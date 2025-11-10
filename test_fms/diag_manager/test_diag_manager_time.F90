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
!! test_diag_manager
!! 2 1 1 0 0 0
!! "ocn%4yr%2mo%2dy%2hr",      1,  "days", 1, "days", "time", 1, "days", "2 1 1 0 0 0"
!! "test_diag_manager_mod", "sst", "sst", "ocn%4yr%2mo%2dy%2hr",  "all", .true., "none", 2

program test_diag_manager_time

use   mpp_domains_mod
use   diag_manager_mod
use   fms_mod
use  time_manager_mod, only: time_type, set_calendar_type, set_date, NOLEAP, JULIAN, operator(+), set_time, print_time

implicit none

type(time_type)                   :: Time
integer, dimension(2)             :: layout = (/1,1/)
integer :: nlon, nlat, nz
type(domain2d)                    :: Domain
real, dimension(:), allocatable :: x, y, z
integer :: i, j
integer :: is, ie, js, je
real, allocatable :: sst(:,:,:), ice(:,:)
integer :: id_x, id_y, id_z, id_sst, id_ice
logical :: used

call fms_init
call set_calendar_type(JULIAN)
call diag_manager_init

nlon = 20
nlat = 30
nz = 5

call mpp_domains_set_stack_size(17280000)
call mpp_define_domains( (/1,nlon,1,nlat/), layout, Domain, name='test_diag_manager')
call mpp_define_io_domain(Domain, (/1,1/))
call mpp_get_compute_domain(Domain, is, ie, js, je)

! Set up the data
allocate(x(nlon), y(nlat), z(nz))
allocate(sst(is:ie,js:je,1:nz), ice(is:ie,js:je))

do i=1,nlon
  x(i) = i
enddo
do j=1,nlat
  y(j) = j
enddo
do i=1,nz
   z(i) = i
enddo

sst = 666.66
ice = 619.0

! Set up the initial time
Time = set_date(2,1,1,0,0,0)

! Register the diags
id_x  = diag_axis_init('x',  x,  'point_E', 'x', long_name='point_E', Domain2=Domain)
id_y  = diag_axis_init('y',  y,  'point_N', 'y', long_name='point_N', Domain2=Domain)
id_z  = diag_axis_init('z',  z,  'point_Z', 'z', long_name='point_Z')
id_sst = register_diag_field  ('test_diag_manager_mod', 'sst', (/id_x,id_y,id_z/), Time, 'SST', 'K')
id_ice = register_diag_field  ('test_diag_manager_mod', 'ice', (/id_x,id_y/), Time, 'ICE', 'm')

! Increase the time and send data
do i=1,23
Time = set_date(2,1,1,i,0,0)
sst = real(i)
ice = real(i)
if(id_sst > 0) used = send_data(id_sst, sst, Time)
if(id_ice > 0) used = send_data(id_ice, ice, Time)
enddo

call diag_manager_end(Time)
call fms_end

end program test_diag_manager_time
