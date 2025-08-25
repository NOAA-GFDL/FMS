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
!> Include file to hold common routines for the coupler tests
!! Uses the FMS_TEST_BC_TYPE_ macro to test both r4/r8

!! constants used for test
#define BCNUM_ 2
#define FLDNUM_ 2

subroutine set_up_1d_coupler_type(bc_type, data_grid)
   type(coupler_1d_bc_type), intent(inout) :: bc_type    !< Coupler 2d restart types
   integer, dimension(2), intent(in)       :: data_grid  !< Starting and ending indexes of data_domain
   integer :: nfiles, nfields, i, j
   character(len=1) :: field_num, file_num

   !bc_type%isc = data_grid(1); bc_type%iec = data_grid(2)
   !bc_type%isd = data_grid(1); bc_type%ied = data_grid(2)
   bc_type%set = .true.
   bc_type%num_bcs = BCNUM_
   nfiles = bc_type%num_bcs
   allocate(bc_type%FMS_TEST_BC_TYPE_(nfiles))

   do i = 1, nfiles
      write(file_num,'(i1)') i
      bc_type%FMS_TEST_BC_TYPE_(i)%ice_restart_file="default_"//file_num//"_ice_restart_2d.nc"

      bc_type%FMS_TEST_BC_TYPE_(i)%num_fields=FLDNUM_
      nfields = bc_type%FMS_TEST_BC_TYPE_(i)%num_fields
      allocate(bc_type%FMS_TEST_BC_TYPE_(i)%field(nfields))

      do j = 1, nfields
         write(field_num,'(i1)') j
         bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%name="var_"//field_num
         allocate(bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%values(data_grid(1):data_grid(2)))
         bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%values = real(j, kind=TEST_FMS_KIND_)
      end do
   end do

end subroutine

!> @brief  Sets up the coupler_2d_bc_type with the information needed to run the test
subroutine set_up_2d_coupler_type(bc_type, data_grid, appendix, to_read)
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

   bc_type%set = .true.
   if (to_read) then
       bc_type%num_bcs = BCNUM_ + 1
   else
       bc_type%num_bcs = BCNUM_
   endif

   bc_type%isc = data_grid(1); bc_type%iec = data_grid(2)
   bc_type%jsc = data_grid(3); bc_type%jec = data_grid(4)
   bc_type%isd = data_grid(1); bc_type%ied = data_grid(2)
   bc_type%jsd = data_grid(3); bc_type%jed = data_grid(4)

   nfiles = bc_type%num_bcs
   allocate(bc_type%FMS_TEST_BC_TYPE_(nfiles))

   do i = 1, nfiles
      write(file_num,'(i1)') i
      if (i==3) then
          bc_type%FMS_TEST_BC_TYPE_(i)%ice_restart_file="default_"//file_num//"_ice_restart_2d.nc"
      else
          bc_type%FMS_TEST_BC_TYPE_(i)%ice_restart_file=appendix//"_"//file_num//"_ice_restart_2d.nc"
      endif

      bc_type%FMS_TEST_BC_TYPE_(i)%num_fields=FLDNUM_
      nfields = bc_type%FMS_TEST_BC_TYPE_(i)%num_fields
      allocate(bc_type%FMS_TEST_BC_TYPE_(i)%field(nfields))

      do j = 1, nfields
         write(field_num,'(i1)') j
         bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%name="var_"//field_num
         allocate(bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%values(data_grid(1):data_grid(2), data_grid(3):data_grid(4)))
         if (to_read) then
             bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%values = real(999., kind=TEST_FMS_KIND_)
         else
             bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%values = real(j, kind=TEST_FMS_KIND_)
         endif
      end do
   end do

endsubroutine

!> @brief Cleans up the coupler_2d_bc_type
subroutine destroy_2d_coupler_type(bc_type)
   type(coupler_2d_bc_type), intent(inout) :: bc_type !< Coupler 2d restart types

   integer :: i,j !< No description

   do i=1, bc_type%num_bcs
      do j = 1, bc_type%FMS_TEST_BC_TYPE_(i)%num_fields
         deallocate(bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%values)
      end do
      deallocate(bc_type%FMS_TEST_BC_TYPE_(i)%field)
   end do

   deallocate(bc_type%FMS_TEST_BC_TYPE_)

end subroutine

!> @brief Compares the data between two coupler_2d_bc_types
subroutine compare_2d_answers(bc_type_read, bc_type)
   type(coupler_2d_bc_type), intent(inout) :: bc_type_read !< Coupler 2d restart types read in
   type(coupler_2d_bc_type), intent(inout) :: bc_type      !< Coupler 2d restart types

   integer :: i,j !< No description

   do i=1, bc_type%num_bcs
      do j = 1, bc_type%FMS_TEST_BC_TYPE_(i)%num_fields
         if (sum(bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%values) .ne. &
             sum(bc_type_read%FMS_TEST_BC_TYPE_(i)%field(j)%values)) then
             call mpp_error(FATAL, "test_coupler_2d: Answers do not match for: " &
                            & //trim(bc_type%FMS_TEST_BC_TYPE_(i)%ice_restart_file)&
                            & //" var: "//trim(bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%name))
         endif
      end do
   end do

   !< Check the dummy general file
   if (sum(bc_type_read%FMS_TEST_BC_TYPE_(3)%field(1)%values) .ne. &
       sum(bc_type_read%FMS_TEST_BC_TYPE_(1)%field(1)%values)) then
       call mpp_error(FATAL, "test_coupler_2d: Answers do not match for var: "// &
                              trim(bc_type_read%FMS_TEST_BC_TYPE_(3)%field(1)%name))
   endif

   if (sum(bc_type_read%FMS_TEST_BC_TYPE_(3)%field(2)%values) .ne. &
       sum(bc_type_read%FMS_TEST_BC_TYPE_(1)%field(2)%values)) then
       call mpp_error(FATAL, "test_coupler_2d: Answers do not match for var: "// &
                              trim(bc_type_read%FMS_TEST_BC_TYPE_(3)%field(2)%name))
   endif


end subroutine

!! 3d routines

!> @brief  Sets up the coupler_3d_bc_type with the information needed to run the test
subroutine set_up_3d_coupler_type(bc_type, data_grid, appendix, to_read)
   type(coupler_3d_bc_type), intent(inout) :: bc_type    !< Coupler 3d restart types
   integer, dimension(5), intent(in)       :: data_grid  !< Starting and ending indexes of data_domain
   character(len=3), intent(in)            :: appendix   !< Appendix to be added to the filename
   logical, intent(in)                     :: to_read    !< If true, it is reading a file so data is
                                                         !! set do a dummy value

   character(len=1) :: field_num  !< string with field_num
   character(len=1) :: file_num   !< string with file_num

   integer :: nfields !< Number of variables
   integer :: nfiles  !< Number of files
   integer :: i,j     !< No description

   bc_type%isc = data_grid(1); bc_type%iec = data_grid(2)
   bc_type%jsc = data_grid(3); bc_type%jec = data_grid(4)
   bc_type%isd = data_grid(1); bc_type%ied = data_grid(2)
   bc_type%jsd = data_grid(3); bc_type%jed = data_grid(4)

   bc_type%set = .true.
   if (to_read) then
       bc_type%num_bcs = BCNUM_ + 1
   else
       bc_type%num_bcs = BCNUM_
   endif

   nfiles = bc_type%num_bcs
   allocate(bc_type%FMS_TEST_BC_TYPE_(nfiles))
   bc_type%ks = 1
   bc_type%ke = data_grid(5)

   do i = 1, nfiles
      write(file_num,'(i1)') i
      if (i==3) then
          bc_type%FMS_TEST_BC_TYPE_(i)%ice_restart_file="default_"//file_num//"_ice_restart_3d.nc"
      else
          bc_type%FMS_TEST_BC_TYPE_(i)%ice_restart_file=appendix//"_"//file_num//"_ice_restart_3d.nc"
      endif

      bc_type%FMS_TEST_BC_TYPE_(i)%num_fields=FLDNUM_
      nfields = bc_type%FMS_TEST_BC_TYPE_(i)%num_fields
      allocate(bc_type%FMS_TEST_BC_TYPE_(i)%field(nfields))

      do j = 1, nfields
         write(field_num,'(i1)') j
         bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%name="var_"//field_num
         allocate(bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%values(data_grid(1):data_grid(2), data_grid(3):data_grid(4), &
                                                               data_grid(5)))

         if (to_read) then
             bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%values = real(999., kind=TEST_FMS_KIND_)
         else
             bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%values = real(j, kind=TEST_FMS_KIND_)
         endif
      end do
   end do

endsubroutine

!> @brief Cleans up the coupler_3d_bc_type
subroutine destroy_3d_coupler_type(bc_type)
   type(coupler_3d_bc_type), intent(inout) :: bc_type !< Coupler 3d restart types

   integer :: i,j !< No description

   do i=1, bc_type%num_bcs
      do j = 1, bc_type%FMS_TEST_BC_TYPE_(i)%num_fields
         deallocate(bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%values)
      end do
      deallocate(bc_type%FMS_TEST_BC_TYPE_(i)%field)
   end do

   deallocate(bc_type%FMS_TEST_BC_TYPE_)

end subroutine

!> @brief Compares the data between two coupler_3d_bc_types
subroutine compare_3d_answers(bc_type_read, bc_type)
   type(coupler_3d_bc_type), intent(inout) :: bc_type_read !< Coupler 3d restart types read in
   type(coupler_3d_bc_type), intent(inout) :: bc_type      !< Coupler 3d restart types

   integer :: i,j !< No description

   do i=1, bc_type%num_bcs
      do j = 1, bc_type%FMS_TEST_BC_TYPE_(i)%num_fields
         if (sum(bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%values) .ne. &
             sum(bc_type_read%FMS_TEST_BC_TYPE_(i)%field(j)%values)) then
             call mpp_error(FATAL, "test_coupler_3d: Answers do not match for: "// &
                            & trim(bc_type%FMS_TEST_BC_TYPE_(i)%ice_restart_file) &
                            & //" var: "//trim(bc_type%FMS_TEST_BC_TYPE_(i)%field(j)%name))
         endif
      end do
   end do

   !< Check the dummy general file
   if (sum(bc_type_read%FMS_TEST_BC_TYPE_(3)%field(1)%values) .ne. &
       sum(bc_type_read%FMS_TEST_BC_TYPE_(1)%field(1)%values)) then
       call mpp_error(FATAL, "test_coupler_3d: Answers do not match for var: "// &
                      & trim(bc_type_read%FMS_TEST_BC_TYPE_(3)%field(1)%name))
   endif

   if (sum(bc_type_read%FMS_TEST_BC_TYPE_(3)%field(2)%values) .ne. &
       sum(bc_type_read%FMS_TEST_BC_TYPE_(1)%field(2)%values)) then
       call mpp_error(FATAL, "test_coupler_3d: Answers do not match for var: "// &
                      & trim(bc_type_read%FMS_TEST_BC_TYPE_(3)%field(2)%name))
   endif


end subroutine
