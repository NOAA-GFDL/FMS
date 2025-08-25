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

!> @brief  This programs tests the functionality in
!! 1. coupler_type_register_restarts (CT_register_restarts_3d)
!! 2. coupler_type_restore_state (CT_restore_state_3d)
program test_coupler_3d

use   fms2_io_mod,        only: FmsNetcdfDomainFile_t, open_file, close_file, read_restart, write_restart
use   fms2_io_mod,        only: FmsNetcdfFile_t, register_axis, register_field, write_data
use   fms2_io_mod,        only: register_variable_attribute
use   fms_mod,            only: fms_init, fms_end
use   mpp_mod,            only: mpp_error, mpp_pe, mpp_root_pe, FATAL
use   mpp_domains_mod,    only: domain2d, mpp_define_domains, mpp_define_io_domain, mpp_get_data_domain, &
                                & mpp_domains_set_stack_size
use   coupler_types_mod,  only: coupler_3d_bc_type, coupler_2d_bc_type, coupler_type_register_restarts, &
                                coupler_type_restore_state
use   coupler_types_mod,  only: coupler_1d_bc_type
use   platform_mod,       only: TEST_FMS_KIND_

implicit none

type(coupler_3d_bc_type)              :: bc_type          !< Coupler 3d restart types
type(coupler_3d_bc_type)              :: bc_type_read     !< Coupler 3d restart types for reading
type(FmsNetcdfDomainFile_t), pointer  :: bc_rest_files(:)=> null() !< Array of fms2_io fileobjs
type(domain2d)                        :: Domain           !< Domain with mask table
integer, dimension(2)                 :: layout = (/1,1/) !< Domain layout
integer                               :: nlon             !< Number of points in x axis
integer                               :: nlat             !< Number of points in y axis
integer, dimension(5)                 :: data_grid        !< Starting/Ending indices in x and y, size of z dimension
                                                          !! for the data_domain
integer                               :: num_rest_files   !< Number of restart files
integer                               :: i                !< No description
type(FmsNetcdfFile_t)                 :: fileobj          !< fms2_io fileobjs
real(TEST_FMS_KIND_), allocatable                     :: dummy_var(:,:,:) !< Dummy variable

call fms_init()

nlat=60
nlon=60

call mpp_domains_set_stack_size(   72000)

call mpp_define_domains( (/1,nlon,1,nlat/), layout, Domain, name='test_coupler')
call mpp_define_io_domain(Domain, (/1,1/))
call mpp_get_data_domain(Domain, data_grid(1), data_grid(2), data_grid(3), data_grid(4))

data_grid(5) = 10

!> Create a dummy general file
if (mpp_pe() .eq. mpp_root_pe()) then
   if (open_file(fileobj, "RESTART/default_3_ice_restart_3d.res.nc", "overwrite")) then
       call register_axis(fileobj, "lonx", nlon)
       call register_axis(fileobj, "laty", nlat)
       call register_axis(fileobj, "z", data_grid(5))

       call register_field(fileobj, "lonx", "double", (/ "lonx" /))
       call register_field(fileobj, "laty", "double", (/ "laty" /))

       call register_field(fileobj, "var_1", "double", (/ "lonx", "laty", "z   " /))
       call register_field(fileobj, "var_2", "double", (/ "lonx", "laty", "z   " /))

       call register_variable_attribute(fileobj, "lonx", "axis", "x", str_len=1)
       call register_variable_attribute(fileobj, "laty", "axis", "y", str_len=1)

       allocate(dummy_var(nlon, nlat, data_grid(5)))
       dummy_var = real(1, kind=TEST_FMS_KIND_)
       call write_data(fileobj, "var_1", dummy_var)

       dummy_var = real(2, kind=TEST_FMS_KIND_)
       call write_data(fileobj, "var_2", dummy_var)

       call close_file(fileobj)

       deallocate(dummy_var)
   endif


endif

!> Write the file with new io
call set_up_3d_coupler_type(bc_type, data_grid, appendix="new", to_read=.false.)
call coupler_type_register_restarts(bc_type, bc_rest_files, num_rest_files, domain, to_read=.false., &
                                  & ocean_restart=.false., directory="RESTART/")

do i = 1, bc_type%num_bcs
   call write_restart(bc_rest_files(i))
   call close_file(bc_rest_files(i))
enddo

!< Now read the file back!
call set_up_3d_coupler_type(bc_type_read, data_grid, appendix="new", to_read=.true.)
call coupler_type_register_restarts(bc_type_read, bc_rest_files, num_rest_files, domain, to_read=.true., &
                                  & ocean_restart=.false., directory="RESTART/")

do i = 1, bc_type_read%num_bcs
   call read_restart(bc_rest_files(i))
enddo

call coupler_type_restore_state(bc_type_read, .true., test_by_field=.true.)

do i = 1, bc_type_read%num_bcs
   call close_file(bc_rest_files(i))
enddo

!< Compare answers!
call compare_3d_answers(bc_type_read, bc_type)

call destroy_3d_coupler_type(bc_type)
call destroy_3d_coupler_type(bc_type_read)

deallocate(bc_rest_files)

call fms_end()

contains

#include "test_coupler_utils.inc"

end program test_coupler_3d
