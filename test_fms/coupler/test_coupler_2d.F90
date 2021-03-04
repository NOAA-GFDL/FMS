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
!! 1. coupler_type_register_restarts (CT_register_restarts_2d)
!! 2. coupler_type_restore_state (CT_restore_state_2d)
program test_coupler_2d

use   fms2_io_mod,        only: FmsNetcdfDomainFile_t, open_file, close_file, read_restart, write_restart
use   fms2_io_mod,        only: FmsNetcdfFile_t, register_axis, register_field, write_data
use   fms2_io_mod,        only: register_variable_attribute
use   fms_mod,            only: fms_init, fms_end
use   mpp_mod,            only: mpp_error, mpp_pe, mpp_root_pe, FATAL
use   mpp_domains_mod,    only: domain2d, mpp_define_domains, mpp_define_io_domain, mpp_get_data_domain
use   coupler_types_mod,  only: coupler_2d_bc_type, coupler_type_register_restarts, coupler_type_restore_state
use   platform_mod,       only: r8_kind

implicit none

type(coupler_2d_bc_type)              :: bc_type          !< Coupler 2d restart types
type(coupler_2d_bc_type)              :: bc_type_read     !< Coupler 2d restart types for reading
type(FmsNetcdfDomainFile_t), pointer  :: bc_rest_files(:)=> null() !< Array of fms2_io fileobjs
type(domain2d)                        :: Domain           !< Domain with mask table
integer, dimension(2)                 :: layout = (/1,1/) !< Domain layout
integer                               :: nlon             !< Number of points in x axis
integer                               :: nlat             !< Number of points in y axis
integer, dimension(4)                 :: data_grid        !< Starting/Ending indices in x and y
                                                          !! for the data_domain
integer                               :: num_rest_files   !< Number of restart files
integer                               :: i                !< No description
type(FmsNetcdfFile_t)                 :: fileobj          !< fms2_io fileobjs
real, allocatable                     :: dummy_var(:,:)   !< Dummy variable

call fms_init()

nlat=60
nlon=60

call mpp_define_domains( (/1,nlon,1,nlat/), layout, Domain, name='test_coupler')
call mpp_define_io_domain(Domain, (/1,1/))
call mpp_get_data_domain(Domain, data_grid(1), data_grid(2), data_grid(3), data_grid(4))

!> Create a dummy general file
if (mpp_pe() .eq. mpp_root_pe()) then
   if (open_file(fileobj, "RESTART/default_3_ice_restart_2d.res.nc", "overwrite")) then
       call register_axis(fileobj, "lonx", nlon)
       call register_axis(fileobj, "laty", nlat)

       call register_field(fileobj, "lonx", "double", (/ "lonx" /))
       call register_field(fileobj, "laty", "double", (/ "laty" /))

       call register_field(fileobj, "var_1", "double", (/ "lonx", "laty" /))
       call register_field(fileobj, "var_2", "double", (/ "lonx", "laty" /))

       call register_variable_attribute(fileobj, "lonx", "axis", "x", str_len=1)
       call register_variable_attribute(fileobj, "laty", "axis", "y", str_len=1)

       allocate(dummy_var(nlon, nlat))
       dummy_var = real(1, kind=r8_kind)
       call write_data(fileobj, "var_1", dummy_var)

       dummy_var = real(2, kind=r8_kind)
       call write_data(fileobj, "var_2", dummy_var)

       call close_file(fileobj)

       deallocate(dummy_var)
   endif


endif


!> Write the file with new io
call set_up_coupler_type(bc_type, data_grid, appendix="new", to_read=.false.)
call coupler_type_register_restarts(bc_type, bc_rest_files, num_rest_files, domain, to_read=.false., ocean_restart=.false., &
                                    & directory="RESTART/")

do i = 1, bc_type%num_bcs
   call write_restart(bc_rest_files(i))
   call close_file(bc_rest_files(i))
enddo

!< Now read the file back!
call set_up_coupler_type(bc_type_read, data_grid, appendix="new", to_read=.true.)
call coupler_type_register_restarts(bc_type_read, bc_rest_files, num_rest_files, domain, to_read=.true., ocean_restart=.false.,&
                                    & directory="RESTART/")

do i = 1, bc_type_read%num_bcs
   call read_restart(bc_rest_files(i))
enddo

call coupler_type_restore_state(bc_type_read, .true., test_by_field=.true.)

do i = 1, bc_type_read%num_bcs
   call close_file(bc_rest_files(i))
enddo

!< Compare answers!
call compare_answers(bc_type_read, bc_type)

call destroy_coupler_type(bc_type)
call destroy_coupler_type(bc_type_read)

deallocate(bc_rest_files)

call fms_end()

contains

!> @brief  Sets up the coupler_2d_bc_type with the information needed to run the test
subroutine set_up_coupler_type(bc_type, data_grid, appendix, to_read)
   type(coupler_2d_bc_type), intent(inout) :: bc_type    !< Coupler 2d restart types
   integer, dimension(4), intent(in)       :: data_grid  !< Starting and ending indexes of data_domain
   character(len=3), intent(in)            :: appendix   !< Appendix to be added to the filename
   logical, intent(in)                     :: to_read    !< If true, it is reading a file so data is
                                                         !! set do a dummy value

   character(len=1) :: field_num  !< string with field_num
   character(len=1) :: file_num   !< string with file_num

   integer :: nfields !< Number of variables
   integer :: nfiles  !< Number of files
   integer :: i,j     !< No description

   if (to_read) then
       bc_type%num_bcs = 3
   else
       bc_type%num_bcs = 2
   endif

   nfiles = bc_type%num_bcs
   allocate(bc_type%bc(nfiles))

   do i = 1, nfiles
      write(file_num,'(i1)') i
      if (i==3) then
          bc_type%bc(i)%ice_restart_file="default_"//file_num//"_ice_restart_2d.nc"
      else
          bc_type%bc(i)%ice_restart_file=appendix//"_"//file_num//"_ice_restart_2d.nc"
      endif

      bc_type%bc(i)%num_fields=2
      nfields = bc_type%bc(i)%num_fields
      allocate(bc_type%bc(i)%field(nfields))

      do j = 1, nfields
         write(field_num,'(i1)') j
         bc_type%bc(i)%field(j)%name="var_"//field_num
         allocate(bc_type%bc(i)%field(j)%values(data_grid(1):data_grid(2), data_grid(3):data_grid(4)))
         if (to_read) then
             bc_type%bc(i)%field(j)%values = real(999., kind=r8_kind)
         else
             bc_type%bc(i)%field(j)%values = real(j, kind=r8_kind)
         endif
      end do
   end do

endsubroutine

!> @brief Cleans up the coupler_2d_bc_type
subroutine destroy_coupler_type(bc_type)
   type(coupler_2d_bc_type), intent(inout) :: bc_type !< Coupler 2d restart types

   integer :: i,j !< No description

   do i=1, bc_type%num_bcs
      do j = 1, bc_type%bc(i)%num_fields
         deallocate(bc_type%bc(i)%field(j)%values)
      end do
      deallocate(bc_type%bc(i)%field)
   end do

   deallocate(bc_type%bc)

end subroutine

!> @brief Compares the data between two coupler_2d_bc_types
subroutine compare_answers(bc_type_read, bc_type)
   type(coupler_2d_bc_type), intent(inout) :: bc_type_read !< Coupler 2d restart types read in
   type(coupler_2d_bc_type), intent(inout) :: bc_type      !< Coupler 2d restart types

   integer :: i,j !< No description

   do i=1, bc_type%num_bcs
      do j = 1, bc_type%bc(i)%num_fields
         if (sum(bc_type%bc(i)%field(j)%values) .ne. sum(bc_type_read%bc(i)%field(j)%values)) then
             call mpp_error(FATAL, "test_coupler_2d: Answers do not match for: "//trim(bc_type%bc(i)%ice_restart_file)//&
                            &" var: "//trim(bc_type%bc(i)%field(j)%name))
         endif
      end do
   end do

   !< Check the dummy general file
   if (sum(bc_type_read%bc(3)%field(1)%values) .ne. sum(bc_type_read%bc(1)%field(1)%values)) then
       call mpp_error(FATAL, "test_coupler_2d: Answers do not match for var: "//trim(bc_type_read%bc(3)%field(1)%name))
   endif

   if (sum(bc_type_read%bc(3)%field(2)%values) .ne. sum(bc_type_read%bc(1)%field(2)%values)) then
       call mpp_error(FATAL, "test_coupler_2d: Answers do not match for var: "//trim(bc_type_read%bc(3)%field(2)%name))
   endif


end subroutine

end program test_coupler_2d
