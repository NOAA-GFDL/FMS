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
! Ryan Mulhall 8/23

!! defaults to ensure compilation 
#ifndef FMS_CP_TEST_KIND_
#define FMS_CP_TEST_KIND_ r8_kind
#endif

#ifndef FMS_TEST_BC_TYPE_
#define FMS_TEST_BC_TYPE_ bc 
#endif

!> Tests for the coupler types interfaces not tested in test_coupler_2d/3d
program test_coupler_types

use   fms_mod,            only: fms_init, fms_end, stdout
use   mpp_mod,            only: mpp_error, mpp_pe, mpp_root_pe, FATAL, mpp_sync
use   mpp_domains_mod,    only: domain2d, mpp_define_domains, mpp_define_io_domain, mpp_get_data_domain, domain1D
use   mpp_domains_mod,    only: mpp_domains_set_stack_size
use   coupler_types_mod,  only: coupler_3d_bc_type, coupler_2d_bc_type, coupler_1d_bc_type
use   coupler_types_mod,  only: coupler_type_copy, coupler_type_spawn, coupler_type_copy_data, coupler_type_redistribute_data
use   coupler_types_mod,  only: coupler_type_rescale_data, coupler_type_increment_data, coupler_type_extract_data, coupler_type_set_data
use   coupler_types_mod,  only: coupler_type_set_diags, coupler_type_write_chksums, coupler_type_send_data, coupler_type_data_override
use   coupler_types_mod,  only: coupler_type_destructor, coupler_type_initialized
use   diag_manager_mod,   only: diag_axis_init, diag_manager_end, diag_manager_init, NULL_AXIS_ID
use   time_manager_mod,   only: time_type, set_date, time_manager_init, set_calendar_type, GREGORIAN
use   platform_mod,       only: r8_kind, r4_kind
implicit none

type(coupler_1d_bc_type) :: bc_1d_new
type(coupler_2d_bc_type) :: bc_2d_new, bc_2d_cp
type(coupler_3d_bc_type) :: bc_3d_new, bc_3d_cp
type(domain2D) :: Domain, Domain_out
integer :: layout(2)
integer :: nlat, nlon, nz, i
integer :: data_grid(5) !< i/j starting and ending indices for data domain
character(len=3) :: appendix !< appoendix added to filename
type(time_type) :: time_t
integer, parameter :: lkind = FMS_CP_TEST_KIND_
real(FMS_CP_TEST_KIND_), allocatable :: array_2d(:,:), array_3d(:,:,:)
integer, parameter :: num_bc = 2, num_fields = 1 !< these are set in set_up_coupler_type routines
real(FMS_CP_TEST_KIND_), allocatable :: lats(:), lons(:), nzs(:) !< arrays of coordinate values for diag_axis initalization
integer :: id_x, id_y, id_z, chksum_unit 
character(len=128) :: chksum_2d, chksum_3d 

call fms_init
call time_manager_init
call set_calendar_type(GREGORIAN)

! basic domain set up
nlat=60; nlon=60; nz=12
layout = (/2, 2/)
call mpp_domains_set_stack_size(86400)
call mpp_define_domains( (/1,nlon,1,nlat/), layout, Domain, name='test_coupler')
call mpp_define_io_domain(Domain, (/1,1/))
call mpp_get_data_domain(Domain, data_grid(1), data_grid(2), data_grid(3), data_grid(4))

! create/allocate new types with routines in utils file
call set_up_1d_coupler_type(bc_1d_new, data_grid)
call set_up_2d_coupler_type(bc_2d_new, data_grid, appendix="new", to_read=.false.)
data_grid(5) = nz 
call set_up_3d_coupler_type(bc_3d_new, data_grid, appendix="new", to_read=.false.)

! test coupler_type_copy, coupler_type_copy_data and coupler_type_destructor
time_t = set_date(1, 1, 1) ! year 0, month 1, day 1
! 1d -> 2d, 3d
call coupler_type_copy(bc_1d_new, bc_2d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), " ", (/ NULL_AXIS_ID /), time_t ) 
call coupler_type_copy(bc_1d_new, bc_3d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), data_grid(5), " ", (/ NULL_AXIS_ID /), time_t ) 
call coupler_type_destructor(bc_2d_cp)
call coupler_type_destructor(bc_3d_cp)
! 2d -> 2d, 3d
call coupler_type_copy(bc_2d_new, bc_2d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), " ", (/ NULL_AXIS_ID /), time_t ) 
call coupler_type_copy(bc_2d_new, bc_3d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), data_grid(5), " ", (/ NULL_AXIS_ID /), time_t ) 
call coupler_type_copy_data(bc_2d_new, bc_2d_cp)
call coupler_type_copy_data(bc_2d_new, bc_3d_cp)
call coupler_type_destructor(bc_2d_cp)
call coupler_type_destructor(bc_3d_cp)
! 3d -> 2d, 3d
call coupler_type_copy(bc_3d_new, bc_2d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), " ", (/ NULL_AXIS_ID /), time_t ) 
call coupler_type_copy(bc_3d_new, bc_3d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), data_grid(5), " ", (/ NULL_AXIS_ID /), time_t ) 
call coupler_type_copy_data(bc_3d_new, bc_3d_cp)
call coupler_type_destructor(bc_2d_cp)
call coupler_type_destructor(bc_3d_cp)

! coupler_type_increment_data
! creates copies to increment into original
call coupler_type_copy(bc_2d_new, bc_2d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), " ", (/ 0 /), time_t ) 
call coupler_type_copy_data(bc_2d_new, bc_2d_cp)
call coupler_type_copy(bc_3d_new, bc_3d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), data_grid(5), " ", (/ 0 /), time_t ) 
call coupler_type_copy_data(bc_3d_new, bc_3d_cp)
call coupler_type_increment_data(bc_2d_new, bc_2d_cp)
call coupler_type_increment_data(bc_3d_new, bc_3d_cp)
call coupler_type_destructor(bc_2d_cp)
call coupler_type_destructor(bc_3d_cp)

! coupler_type_rescale_data
call coupler_type_rescale_data(bc_2d_new, 2.0_lkind) 
call coupler_type_rescale_data(bc_3d_new, 2.0_lkind)

! coupler_type_extract_data
! loop through set number of boundary conditions
! using 1 for field_index(?) 
allocate(array_2d(data_grid(1):data_grid(2), data_grid(3):data_grid(4)))
allocate(array_3d(data_grid(1):data_grid(2), data_grid(3):data_grid(4), data_grid(5)))
print *, SIZE(array_2d,1), SIZE(array_2d,2)
do i=1, num_bc
  call coupler_type_extract_data(bc_2d_new, i, num_fields, array_2d)
  call coupler_type_extract_data(bc_3d_new, i, num_fields, array_3d)
enddo

! coupler_type_set_data
array_2d = 1.0_lkind
array_3d = 1.0_lkind
do i=1, num_bc
  call coupler_type_set_data(array_2d, i, num_fields, bc_2d_new)
  call coupler_type_set_data(array_2d, i, num_fields, data_grid(5), bc_3d_new)
  call coupler_type_set_data(array_3d, i, num_fields, bc_3d_new)
enddo



! diag_manager wrapper routines
! coupler_type_set_diags and coupler_type_send_data

! set up for diag manager
call diag_manager_init
allocate(lats(1:nlat), lons(1:nlon), nzs(1:nz))
do i=1, nlat
  lats(i) = i
enddo
do i=1, nlon
  lons(i) = i
enddo
do i=1, nz 
  nzs(i) = i
enddo
id_x  = diag_axis_init('x',  lats,  'point_E', 'x', long_name='point_E', Domain2=Domain)
id_y  = diag_axis_init('y',  lons,  'point_N', 'y', long_name='point_N', Domain2=Domain)
id_z  = diag_axis_init('z',  nzs,  'point_Z', 'z', long_name='point_Z')

! registers field with type
call coupler_type_set_diags(bc_2d_new, "dat1", (/id_x, id_y/), time_t)
call coupler_type_set_diags(bc_3d_new, "dat2", (/id_x, id_y, id_z/), time_t)
! send some data to the fields
call coupler_type_copy(bc_2d_new, bc_2d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), &
                       "dat1", (/ id_x, id_y /), time_t ) 
call coupler_type_copy_data(bc_2d_new, bc_2d_cp)
call coupler_type_copy(bc_3d_new, bc_3d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), data_grid(5), &
                       "dat2", (/ id_x, id_y, id_z /), time_t ) 
call coupler_type_copy_data(bc_3d_new, bc_3d_cp)
do i=1,12
  time_t = set_date(1, i, 1)
  call coupler_type_increment_data(bc_2d_cp, bc_2d_new)
  call coupler_type_increment_data(bc_3d_cp, bc_3d_new)
  call coupler_type_send_data(bc_2d_new, time_t) 
  call coupler_type_send_data(bc_3d_new, time_t) 
enddo
time_t = set_date(2, 1, 1)
call diag_manager_end(time_t)

! coupler_type_write_chksums
! writes to a new ascii file, easier to check than std*'s 
open(newunit=chksum_unit, action='WRITE', file="coupler_chksum.out" )
call coupler_type_write_chksums(bc_2d_new, chksum_unit)
call coupler_type_write_chksums(bc_3d_new, chksum_unit)
close(chksum_unit)
open(newunit=chksum_unit, action='READ', file="coupler_chksum.out" )
read(chksum_unit, '(A)') chksum_2d
print *, '2D checksums'
print *, chksum_2d
call mpp_sync()
print *, '3D checksums'
read(chksum_unit, '(A)') chksum_3d
print *, chksum_3d
close(chksum_unit)

! coupler_type_data_override 
! TODO
call diag_manager_init(Ocean_domain_in=Domain)


! coupler_type_redistribute_data
! just using the same domain
call mpp_define_domains((/1, nlon, 1, nlat/), layout, Domain_out, name="test_coupler_redistributed_2x2")
call set_up_2d_coupler_type(bc_2d_cp, data_grid, appendix="new", to_read=.false.)
call set_up_3d_coupler_type(bc_3d_cp, data_grid, appendix="new", to_read=.false.)
call coupler_type_redistribute_data(bc_2d_new, Domain, bc_2d_cp, domain_out, complete=.true.)
call coupler_type_redistribute_data(bc_3d_new, Domain, bc_3d_cp, domain_out, complete=.true.)
call coupler_type_destructor(bc_2d_cp)
call coupler_type_destructor(bc_3d_cp)
! using a different layout
call mpp_define_domains((/1, nlon, 1, nlat/), (/1, 4/), Domain_out, name="test_coupler_redistributed_1x4")
call mpp_get_data_domain(Domain_out, data_grid(1), data_grid(2), data_grid(3), data_grid(4))
call set_up_2d_coupler_type(bc_2d_cp, data_grid, appendix="new", to_read=.false.)
call set_up_3d_coupler_type(bc_3d_cp, data_grid, appendix="new", to_read=.false.)
call coupler_type_redistribute_data(bc_2d_new, Domain, bc_2d_cp, domain_out, complete=.true.)
call coupler_type_redistribute_data(bc_3d_new, Domain, bc_3d_cp, domain_out, complete=.true.)

! clean up
call coupler_type_destructor(bc_1d_new)
call coupler_type_destructor(bc_2d_new)
call coupler_type_destructor(bc_3d_new)

call fms_end

contains

#include "test_coupler_utils.inc"

end program