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
!----------
!ug support

!>Return a flag indicating whether the inputted field exists in the inputted
!!file, where the file is associated with an unstructured mpp domain.
function fms_io_unstructured_field_exist(file_name, &
                                         field_name, &
                                         domain) &
                                         result(does_field_exist)

   !Inputs/Outputs
    character(len=*),intent(in) :: file_name        !<Name of a file.
    character(len=*),intent(in) :: field_name       !<Name of a field.
    type(domainUG),intent(in)   :: domain           !<An unstructured mpp domain.
    logical(INT_KIND)           :: does_field_exist !<Flag telling if the inputted field exists in the inputted file.

   !Local variables
    logical(INT_KIND)                        :: file_exist !<Flag telling if the inputted file or one of its variants exists.
    character(len=256)                       :: fname      !<Actual name of the found file.
    logical(INT_KIND)                        :: read_dist  !<Flag telling if the file is "distributed" (has IO domain tile id appended onto the end).
    integer(INT_KIND)                        :: funit      !<A file unit.
    integer(INT_KIND)                        :: nfile      !<Index of the inputted file in the "files_read" module array.
    integer(INT_KIND)                        :: i          !<Loop variable.
    integer(INT_KIND)                        :: ndim       !<Number of dimensions in a file.
    integer(INT_KIND)                        :: nvar       !<Number of fields in a file.
    integer(INT_KIND)                        :: natt       !<Number of attributes in a file.
    integer(INT_KIND)                        :: ntime      !<Number of time levels in a file.
    character(len=64)                        :: tmp_name   !<Name of a field.
    type(fieldtype),dimension(:),allocatable :: fields     !<An array of fields found in the input file.


   !Set default return value for the function.
    does_field_exist = .false.

   !Return if the inputted field name is in valid.
    if (len_trim(field_name) .eq. 0) then
        return
    endif
    if (field_name(1:1) .eq. ' ') then
        return
    endif

   !Check if the file exists.
    file_exist = fms_io_unstructured_get_file_name(file_name, &
                                                   fname, &
                                                   read_dist, &
                                                   domain)

    if (file_exist) then

       !Get the file unit for the input file.
        call fms_io_unstructured_get_file_unit(fname, &
                                               funit, &
                                               nfile, &
                                               read_dist, &
                                               domain)

       !Get the number of dimensions, fields, attributes and time levels
       !for the file.
        call mpp_get_info(funit, &
                          ndim, &
                          nvar, &
                          natt, &
                          ntime)

       !Create an array of all fields contained in the file.
        allocate(fields(nvar))
        call mpp_get_fields(funit, &
                            fields)

       !Loop through the fields to see if the inputted field name matches
       !any of the fields from the file.
        do i = 1,nvar
            call mpp_get_atts(fields(i), &
                              name=tmp_name)
            if (lowercase(trim(tmp_name)) .eq. lowercase(trim(field_name))) then
                does_field_exist = .true.
            endif
        enddo

       !Deallocate local allocatable.
        deallocate(fields)
    endif

    return
end function fms_io_unstructured_field_exist

!----------
