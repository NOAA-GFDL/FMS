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

use   fms_mod,            only: fms_init, fms_end, stdout, string
use   mpp_mod,            only: mpp_error, mpp_pe, mpp_root_pe, FATAL, mpp_sync, mpp_init, input_nml_file
use   mpp_domains_mod,    only: domain2d, mpp_define_domains, mpp_define_io_domain, mpp_get_data_domain, domain1D
use   mpp_domains_mod,    only: mpp_domains_set_stack_size
use   coupler_types_mod,  only: coupler_3d_bc_type, coupler_2d_bc_type, coupler_1d_bc_type
use   coupler_types_mod,  only: coupler_type_copy, coupler_type_spawn, coupler_type_copy_data
use   coupler_types_mod,  only: coupler_type_redistribute_data, coupler_type_set_data, coupler_type_data_override
use   coupler_types_mod,  only: coupler_type_rescale_data, coupler_type_increment_data, coupler_type_extract_data
use   coupler_types_mod,  only: coupler_type_set_diags, coupler_type_write_chksums, coupler_type_send_data
use   coupler_types_mod,  only: coupler_type_destructor, coupler_type_initialized
use   diag_manager_mod,   only: diag_axis_init, diag_manager_end, diag_manager_init, NULL_AXIS_ID
use   time_manager_mod,   only: time_type, set_date, time_manager_init, set_calendar_type, JULIAN
use   data_override_mod,  only: data_override_init
use   constants_mod,      only: pi
use   platform_mod,       only: r8_kind, r4_kind
use   fms2_io_mod,        only: fms2_io_init
use   netcdf,             only: nf90_close, nf90_put_var, nf90_enddef, nf90_create, nf90_def_dim, nf90_clobber, &
                                nf90_64bit_offset, nf90_char, nf90_def_var, nf90_float
implicit none

type(coupler_1d_bc_type) :: bc_1d_new
type(coupler_2d_bc_type) :: bc_2d_new, bc_2d_cp
type(coupler_3d_bc_type) :: bc_3d_new, bc_3d_cp
type(coupler_2d_bc_type) :: bc_2d_ref !< just used to check answers
type(coupler_3d_bc_type) :: bc_3d_ref !< just used to check answers
type(domain2D) :: Domain, Domain_out
integer :: layout(2)
integer :: nlat, nlon, nz, i, j
integer :: data_grid(5) !< i/j starting and ending indices for data domain
character(len=3) :: appendix !< appoendix added to filename
type(time_type) :: time_t
integer, parameter :: lkind = FMS_CP_TEST_KIND_
real(FMS_CP_TEST_KIND_), allocatable :: array_2d(:,:), array_3d(:,:,:)
integer, parameter :: num_bc = 2, num_fields = 2 !< these are set in set_up_coupler_type routines
real(FMS_CP_TEST_KIND_), allocatable :: lats(:), lons(:), nzs(:) !< arrays of coordinate values for diag_axis
                                                                 !! initalization
integer :: id_x, id_y, id_z, chksum_unit
character(len=128) :: chksum_2d, chksum_3d
real(FMS_CP_TEST_KIND_), allocatable :: expected_2d(:,:), expected_3d(:,:,:)
integer :: err, ncid, dim1D, varid, day
logical, allocatable :: return_stats(:,:)

logical :: fail_return_status = .false. !< if true checks for one of the coupler_type_send_data calls to fail and
                                        !! return a false value

NAMELIST /test_coupler_types_nml/ fail_return_status

call fms_init
call time_manager_init
call fms2_io_init
call mpp_init
call set_calendar_type(JULIAN)

read(input_nml_file, test_coupler_types_nml, iostat=err)
if(err > 0) call mpp_error(FATAL, "test_coupler_types:: error reading test input nml")

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

! coupler_type_set_data
allocate(array_2d(data_grid(1):data_grid(2), data_grid(3):data_grid(4)))
allocate(array_3d(data_grid(1):data_grid(2), data_grid(3):data_grid(4), data_grid(5)))
array_2d = 1.0_lkind
array_3d = 1.0_lkind
do i=1, num_bc
  do j=1, num_fields
    call coupler_type_set_data(array_2d, i, j, bc_2d_new)
    call coupler_type_set_data(array_2d, i, j, data_grid(5), bc_3d_new)
    call coupler_type_set_data(array_3d, i, j, bc_3d_new)
  enddo
enddo
call check_field_data_2d(bc_2d_new, array_2d)
call check_field_data_3d(bc_3d_new, array_3d)

! coupler_type_write_chksum
! needs to write to a unit num
call coupler_type_write_chksums(bc_2d_new, stdout())
call coupler_type_write_chksums(bc_3d_new, stdout())

! coupler_type_increment_data
! creates copies to increment into original
call coupler_type_copy(bc_2d_new, bc_2d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), " ", &
                       (/ 0 /), time_t )
call coupler_type_copy_data(bc_2d_new, bc_2d_cp)
call coupler_type_copy(bc_3d_new, bc_3d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), data_grid(5), " ", &
                       (/ 0 /), time_t )
call coupler_type_copy_data(bc_3d_new, bc_3d_cp)
call coupler_type_increment_data(bc_2d_new, bc_2d_cp)
call coupler_type_increment_data(bc_3d_new, bc_3d_cp)
! copy of itself incremented should just be 2.0
array_2d = 2.0_lkind; array_3d = 2.0_lkind
call check_field_data_2d(bc_2d_cp, array_2d)
call check_field_data_3d(bc_3d_cp, array_3d)

! coupler_type_rescale_data
call coupler_type_rescale_data(bc_2d_cp, 2.0_lkind)
call coupler_type_rescale_data(bc_3d_cp, 2.0_lkind)
array_2d = 4.0_lkind; array_3d = 4.0_lkind ! data was 2, rescaled by factor of 2
call check_field_data_2d(bc_2d_cp, array_2d)
call check_field_data_3d(bc_3d_cp, array_3d)
call coupler_type_destructor(bc_2d_cp)
call coupler_type_destructor(bc_3d_cp)

! coupler_type_extract_data
do i=1, num_bc
  do j=1, num_fields
    call coupler_type_extract_data(bc_2d_new, i, j, array_2d)
    call coupler_type_extract_data(bc_3d_new, i, j, array_3d)
  enddo
enddo
call check_field_data_2d(bc_2d_new, array_2d)
call check_field_data_3d(bc_3d_new, array_3d)

! test coupler_type_copy, coupler_type_copy_data and coupler_type_destructor
time_t = set_date(1, 1, 1)
! 1d -> 2d, 3d
call coupler_type_copy(bc_1d_new, bc_2d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), " ", &
                       (/ NULL_AXIS_ID /), time_t )
call coupler_type_copy(bc_1d_new, bc_3d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), data_grid(5)," ",&
                       (/ NULL_AXIS_ID /), time_t )
call coupler_type_destructor(bc_2d_cp)
call coupler_type_destructor(bc_3d_cp)
! 2d -> 2d, 3d
call coupler_type_copy(bc_2d_new, bc_2d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), " ", &
                       (/ NULL_AXIS_ID /), time_t )
call coupler_type_copy(bc_2d_new, bc_3d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), data_grid(5), " ", &
                       (/ NULL_AXIS_ID /), time_t )
call coupler_type_copy_data(bc_2d_new, bc_2d_cp)
call coupler_type_copy_data(bc_2d_new, bc_3d_cp)
array_2d = 1.0; array_3d = 1.0
call check_field_data_2d(bc_2d_cp, array_2d)
call check_field_data_3d(bc_3d_cp, array_3d)
call coupler_type_destructor(bc_2d_cp)
call coupler_type_destructor(bc_3d_cp)
! 3d -> 2d, 3d
call coupler_type_copy(bc_3d_new, bc_2d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), " ", &
                       (/ NULL_AXIS_ID /), time_t )
call coupler_type_copy(bc_3d_new, bc_3d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), data_grid(5), " ", &
                       (/ NULL_AXIS_ID /), time_t )
call coupler_type_copy_data(bc_3d_new, bc_3d_cp)
call check_field_data_3d(bc_3d_cp, array_3d)
call coupler_type_destructor(bc_2d_cp)
call coupler_type_destructor(bc_3d_cp)

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
! registers field with data in type
! reset the time and assign names to each field
time_t = set_date(1, 1, 1)
do i=1, num_bc
  do j=1, num_fields
    bc_2d_new%FMS_TEST_BC_TYPE_(i)%field(j)%name = "bc"//string(i)//"_var2d_"//string(j)
    bc_3d_new%FMS_TEST_BC_TYPE_(i)%field(j)%name = "bc"//string(i)//"_var3d_"//string(j)
    bc_2d_new%FMS_TEST_BC_TYPE_(i)%field(j)%long_name = "bc"//string(i)//"_variable_2d_"//string(j)//"_min"
    bc_3d_new%FMS_TEST_BC_TYPE_(i)%field(j)%long_name = "bc"//string(i)//"_variable_3d_"//string(j)//"_min"
  enddo
enddo
call coupler_type_set_diags(bc_2d_new, "test_coupler_types", (/id_x, id_y/), time_t)
call coupler_type_set_diags(bc_3d_new, "test_coupler_types", (/id_x, id_y, id_z/), time_t)
call coupler_type_copy(bc_2d_new, bc_2d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), " ", &
                       (/null_axis_id/), time_t)
call coupler_type_copy_data(bc_2d_new, bc_2d_cp)
call coupler_type_copy(bc_3d_new, bc_3d_cp, data_grid(1), data_grid(2), data_grid(3), data_grid(4), data_grid(5),  " ",&
                       (/null_axis_id/), time_t)
call coupler_type_copy_data(bc_3d_new, bc_3d_cp)

do day=1,31
  time_t = set_date(1, 1, day)
  call coupler_type_increment_data(bc_2d_cp, bc_2d_new) ! increment _new with cp
  call coupler_type_increment_data(bc_3d_cp, bc_3d_new)
  call coupler_type_send_data(bc_2d_new, time_t, return_stats)
  if( fail_return_status ) then
    if( ALL(return_stats) ) call mpp_error(FATAL, "test_coupler_types:: send_data calls returned true, "// &
                                                  "expected false return value from incorrect diag_table")
  else
    if( .not. ALL(return_stats) ) call mpp_error(FATAL, &
                                  "test_coupler_types:: coupler_type_send_data returned false with valid diag_table")
  endif
  call coupler_type_send_data(bc_3d_new, time_t, return_stats)
  if( fail_return_status ) then
    if( ALL(return_stats) ) call mpp_error(FATAL, "test_coupler_types:: send_data calls returned true, "// &
                                                  "expected false return value from incorrect diag_table")
  else
    if( .not. ALL(return_stats) ) call mpp_error(FATAL, &
                                  "test_coupler_types:: coupler_type_send_data returned false with valid diag_table")
  endif
enddo
time_t = set_date(1, 2, 1)
call diag_manager_end(time_t)

! coupler_type_data_override
! basic grid spec points to outputted .nc's
if( mpp_pe() .eq. mpp_root_pe()) then
  err = nf90_create('INPUT/grid_spec.nc', ior(nf90_clobber, nf90_64bit_offset), ncid)
  err = nf90_def_dim(ncid, 'str', 60, dim1d)
  err = nf90_def_var(ncid, 'x_T', nf90_char, (/dim1d/), varid)
  err = nf90_put_var(ncid, varid, "coupler_types_bc1.nc")
  err = nf90_def_var(ncid, 'xta', nf90_float, (/dim1d/), varid)
  err = nf90_def_var(ncid, 'yta', nf90_float, (/dim1d/), varid)
  err = nf90_enddef(ncid)
  err = nf90_close(ncid)
endif
call mpp_sync()
call data_override_init(Atm_domain_in=Domain, mode=FMS_CP_TEST_KIND_)

time_t = set_date(1, 1, 15)
call coupler_type_data_override("ATM", bc_2d_new, time_t)
call coupler_type_data_override("ATM", bc_3d_new, time_t)
call coupler_type_data_override("OCN", bc_2d_new, time_t)
call coupler_type_data_override("OCN", bc_3d_new, time_t)
call coupler_type_data_override("ICE", bc_2d_new, time_t)
call coupler_type_data_override("ICE", bc_3d_new, time_t)
call coupler_type_data_override("LND", bc_2d_new, time_t)
call coupler_type_data_override("LND", bc_3d_new, time_t)

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
! check deallocation
! both should be deallocated regardless of kind
if( associated(bc_1d_new%bc) .or. associated(bc_2d_new%bc) .or. associated(bc_3d_new%bc)) &
  call mpp_error(FATAL, "test_coupler_types: bc type still associated after destructor called")
if( associated(bc_1d_new%bc_r4) .or. associated(bc_2d_new%bc_r4) .or. associated(bc_3d_new%bc_r4)) &
  call mpp_error(FATAL, "test_coupler_types: bc_r4 type still associated after destructor called")

call fms_end

contains

#include "test_coupler_utils.inc"

subroutine check_field_data_2d(bc_2d, expected)
  type(coupler_2d_bc_type) :: bc_2d
  real(FMS_CP_TEST_KIND_), intent(in) :: expected(:,:)
  real(FMS_CP_TEST_KIND_), pointer :: values_ptr(:,:)

  do i=1, bc_2d%num_bcs
    do j=1, bc_2d%FMS_TEST_BC_TYPE_(i)%num_fields
      values_ptr => bc_2d%FMS_TEST_BC_TYPE_(i)%field(j)%values
      ! checks each index
      if(SUM(values_ptr) .ne. SUM(expected)) then
        print *, "SUMS", SUM(values_ptr), SUM(expected), SHAPE(values_ptr), SHAPE(expected)
        call mpp_error(FATAL, "test_coupler_types: incorrect 2d values against expected result")
      endif
    enddo
  enddo
end subroutine

subroutine check_field_data_3d(bc_3d, expected)
  type(coupler_3d_bc_type) :: bc_3d
  real(FMS_CP_TEST_KIND_), intent(in) :: expected(:,:,:)
  real(FMS_CP_TEST_KIND_), pointer :: values_ptr(:,:,:)
  integer :: x, y, z, vals_start(3) !< need start point for indices, passed in will always be 1-n

  do i=1, bc_3d%num_bcs
    do j=1, bc_3d%FMS_TEST_BC_TYPE_(i)%num_fields
      values_ptr => bc_3d%FMS_TEST_BC_TYPE_(i)%field(j)%values
      if(SUM(values_ptr) .ne. SUM(expected)) then
        print *, "SUMS", SUM(values_ptr), SUM(expected), SHAPE(values_ptr), SHAPE(expected)
        call mpp_error(FATAL, "test_coupler_types: incorrect 3d values against expected result")
      endif
    enddo
  enddo
end subroutine check_field_data_3d

end program
