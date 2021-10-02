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

program test_axis_utils

use fms_mod,         only : fms_init, fms_end, check_nml_error
use mpp_mod,         only : mpp_sync, mpp_pe, mpp_root_pe, mpp_error, FATAL, stdout, &
                            mpp_get_current_pelist, mpp_npes
use mpp_mod,         only : input_nml_file
use axis_utils2_mod, only : axis_edges
use fms2_io_mod,     only : open_file, close_file, write_data, register_axis, register_field, &
                            FmsNetcdfFile_t, register_variable_attribute
use platform_mod,    only : r8_kind

implicit none

type data_type
   real(kind=r8_kind) :: var(10)         !< Axis data
   real(kind=r8_kind) :: var_edges(2,10) !< The boundaries of the axis data
   real(kind=r8_kind) :: answers(11)     !< The expected result
end type data_type

type(data_type)       :: data_in     !< Data used to create the netcdf file
integer, allocatable  :: pes(:)      !< List of pes
type(FmsNetcdfFile_t) :: fileobj     !< FMS2_io fileobj

real(kind=r8_kind)    :: answers(11) !< Results obtained from the axis_edges call

call fms_init

!< Get the current pelist
allocate(pes(mpp_npes()))
call mpp_get_current_pelist(pes)

call set_data(data_in)
call create_input_files(data_in)

!< Test calls to axis_edges
if ( .not. open_file(fileobj, "test_axis_utils.nc", "read", pelist=pes)) then
    call mpp_error(FATAL, "Error opening test_axis_utils.nc to read")
endif

!< Case 1: Here the variable "axis" in the file does not have the attribute "bounds" or "edges", so
!! it calculates them from the data in "axis"
answers = 0.0
call axis_edges(fileobj, "axis", answers)
call compare_answers(answers, data_in%answers, "1")

!< Case 2: Here the variable "axis_with_bounds" in the file has the attribute
!! "bounds", so the data is read from the variable "bounds"
answers = 0.0
call axis_edges(fileobj, "axis_with_bounds", answers)
call compare_answers(answers, data_in%answers, "2")

!< Case 3: Here the variable "axis_with_edges" in the file has the attribute
!"edges", so the data is read from the variable "edges"
answers = 0.0
call axis_edges(fileobj, "axis_with_edges", answers)
call compare_answers(answers, data_in%answers, "3")

!< Case 4: Here the flag "reproduce_null_char_bug_flag" is turned on, so the
!! edges are calculated from the data in axis because edges has a null character
!! in the end
answers = 0.0
call axis_edges(fileobj, "axis_with_edges", answers, reproduce_null_char_bug_flag=.true.)
call compare_answers(answers, data_in%answers, "4")

call close_file(fileobj)
deallocate(pes)

call fms_end

contains

!> @brief  Compares the values of two arrays
subroutine compare_answers(answers_in, answers_expected, test_case)
real(kind=r8_kind), intent(in) :: answers_in(:) !< Answer calculated
real(kind=r8_kind), intent(in) :: answers_expected(:) !< Answer expected
character(1),       intent(in) :: test_case !< String indicating the case number

integer :: i !< For do loop

do i = 1, size(answers_expected,1)
   if(answers_in(i) .ne. answers_expected(i)) then
      print *, "i=", i, " Answer in: ", answers_in(i), " Answer expected ", answers_expected(i)
      call mpp_error(FATAL, "axis_edges case"//trim(test_case)//": Answers are not correct")
   endif
enddo
end subroutine compare_answers

!> @brief  Sets the values of the data_type to be use to write the file, and to
!! compare answers
subroutine set_data(data_in)
type(data_type), intent(out) :: data_in !< data_type to set the expected values to

integer :: i !< For do loop

do i=1,10
   data_in%var(i) = real(i, kind=r8_kind)-0.5_r8_kind

   data_in%var_edges(1,i) = real(i-1, kind=r8_kind)
   data_in%var_edges(2,i) = real(i, kind=r8_kind)

   data_in%answers(i) = real(i-1, kind=r8_kind)
enddo

data_in%answers(11) = real(10, kind=r8_kind)

end subroutine

!> @brief  Creates a netcdf file to test the different test cases of
!!"axis_edges"
subroutine create_input_files(data_in)
type(data_type), intent(in) :: data_in !< data_type containing the values to be added to the file

type(FmsNetcdfFile_t) :: fileobj !< FMS2_io fileobj

if (mpp_pe() .eq. mpp_root_pe()) then
   if ( .not. open_file(fileobj, "test_axis_utils.nc", "overwrite")) then
      call mpp_error(FATAL, "Error opening test_axis_utils.nc to write")
   endif

   call register_axis(fileobj, "dim1", 10)
   call register_axis(fileobj, "dim2", 2)

   call register_field(fileobj, "axis", "double", dimensions=(/"dim1"/))

   call register_field(fileobj, "axis_with_bounds", "double", dimensions=(/"dim1"/))
   call register_variable_attribute(fileobj, "axis_with_bounds", "bounds", "bounds", str_len=6)
   call register_field(fileobj, "bounds", "double", dimensions=(/"dim2", "dim1"/))

   call register_field(fileobj, "axis_with_edges", "double", dimensions=(/"dim1"/))
   call register_variable_attribute(fileobj, "axis_with_edges", "edges", "edges"//char(0), str_len=6)
   call register_field(fileobj, "edges", "double", dimensions=(/"dim2", "dim1"/))

   call write_data(fileobj, "axis", data_in%var)
   call write_data(fileobj, "axis_with_bounds", data_in%var)
   call write_data(fileobj, "axis_with_edges", data_in%var)
   call write_data(fileobj, "bounds", data_in%var_edges)
   call write_data(fileobj, "edges", data_in%var_edges)

   call close_file(fileobj)
endif

!< Wait for root_pe to catch up!
call mpp_sync()

end subroutine create_input_files

end program test_axis_utils
