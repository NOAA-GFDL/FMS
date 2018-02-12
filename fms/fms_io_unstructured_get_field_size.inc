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

!>Get the size of the dimensions of a field from a file associated with an
!!unstructured mpp domain.
subroutine fms_io_unstructured_get_field_size(filename, &
                                              fieldname, &
                                              field_dimension_sizes, &
                                              domain, &
                                              field_found)

   !Inputs/Outputs
    character(len=*),intent(in)        :: filename              !<The name of a file.
    character(len=*),intent(in)        :: fieldname             !<The name of a field in the input file.
    integer,dimension(:),intent(inout) :: field_dimension_sizes !<Array of dimension sizes for the inputted field.
    type(domainUG),intent(in)          :: domain                !<An unstructured mpp domain associated with the input file.
    logical,intent(out),optional       :: field_found           !<Flag telling if the inputted field was found in the inputted file.

   !Local variables
    type(domainUG),pointer                     :: io_domain       !<Pointer to the unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes  !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist          !<A pelist.
    integer(INT_KIND)                          :: funit           !<File unit for the inputted file.
    integer(INT_KIND)                          :: num_axes        !<The total number of axes contained in the inputted file.
    integer(INT_KIND)                          :: num_fields      !<The total number of fields contained in the inputted file.
    integer(INT_KIND)                          :: num_atts        !<The total number of global attributes contained in the inputted file.
    integer(INT_KIND)                          :: num_time_levels !<The total number of time levels contained in the inputted file.
    type(fieldtype),dimension(max_fields)      :: file_fields     !<An array of all fields contained in the inputted file (max_fields is a module variable).
    logical(INT_KIND)                          :: found           !<Flag telling if the field was found in the file.
    character(len=128)                         :: file_field_name !<Name of a field from the inputted file.
    integer(INT_KIND)                          :: file_field_ndim !<Number of dimensions of a field from the inputted file.
    type(axistype),dimension(max_fields)       :: file_field_axes !<An array of all axes of a field contained in the inputted file (max_fields is a module variable).
    character(len=128)                         :: file_axis_name  !<Name of an axis from the inputted file.
    integer(INT_KIND)                          :: file_axis_size  !<Size of an axis from the inputted file.
    integer(INT_KIND)                          :: i               !<Loop variable.
    integer(INT_KIND)                          :: j               !<Loop variable.

   !Point to the I/O domain associated with the inputted unstructured mpp
   !domain.
    io_domain => null()
    io_domain => mpp_get_UG_io_domain(domain)

   !Get the pelist associated with the I/O domain.
    io_domain_npes = mpp_get_UG_domain_npes(io_domain)
    allocate(pelist(io_domain_npes))
    call mpp_get_UG_domain_pelist(io_domain, &
                                  pelist)
    io_domain => null()

   !Get the file unit for the inputted file.
    call fms_io_unstructured_file_unit(filename, &
                                       funit, &
                                       domain)

   !Have the root rank of the I/O domain pelist get the size of the dimensions
   !of the inputted fields from the inputted file.
    if (mpp_pe() .eq. pelist(1)) then

       !Get the number of fields and axes contained in the inputted file.
        call mpp_get_info(funit, &
                          num_axes, &
                          num_fields, &
                          num_atts, &
                          num_time_levels)

       !Make sure that the number of fields in the file does not exceed the
       !maximum number allowed per file.
       !max_fields is a module variable.
        if (num_fields .gt. max_fields) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_get_field_size:" &
                           //" the number of fields in the file " &
                           //trim(filename)//" exceeds the maximum number" &
                           //" of fields allowed per file (max_fields)")
        endif

       !Read in all fields contained in the inputted file.
        call mpp_get_fields(funit, &
                            file_fields(1:num_fields))

       !Check if the inputted field matches one the fields contained in
       !the inputted file.  If it matches, get the size of the field
       !dimensions.
        found = .false.
        field_dimension_sizes = -1
        do i = 1,num_fields
            call mpp_get_atts(file_fields(i), &
                              name=file_field_name)
            if (lowercase(trim(file_field_name)) .eq. &
                lowercase(trim(fieldname))) then
                call mpp_get_atts(file_fields(i), &
                                  ndim=file_field_ndim)
                call mpp_get_atts(file_fields(i), &
                                  axes=file_field_axes(1:file_field_ndim))
                do j = 1,file_field_ndim
                    call mpp_get_atts(file_field_axes(j), &
                                      len=field_dimension_sizes(j))
                enddo
                found = .true.
                exit
            endif
        enddo

       !If the inputted field does not match any of the fields contained
       !in the inputted file, then check if it matches any of the axes
       !contained in the file.
        if (.not. found) then
            call mpp_get_axes(funit, &
                              file_field_axes(1:num_axes))
            do i = 1,num_axes
                call mpp_get_atts(file_field_axes(i), &
                                  name=file_axis_name, &
                                  len=file_axis_size)
                if (lowercase(trim(file_axis_name)) .eq. &
                    lowercase(trim(fieldname))) then
                    field_dimension_sizes(1) = file_axis_size
                    found = .true.
                    exit
                endif
            enddo
        endif
    endif

   !Broadcast the flag telling if the inputted field was found in the inputted
   !file and the field dimension sizes array to all non-root ranks on the
   !I/O domain pelist.
    if (mpp_pe() .eq. pelist(1)) then
        do i = 2,io_domain_npes
            call mpp_send(found, &
                          pelist(i), &
                          tag=COMM_TAG_1)
            call mpp_send(field_dimension_sizes, &
                          size(field_dimension_sizes), &
                          pelist(i), &
                          tag=COMM_TAG_2)
        enddo
        call mpp_sync_self()
    else
        call mpp_recv(found, &
                      pelist(1), &
                      block = .false., &
                      tag=COMM_TAG_1)
        call mpp_recv(field_dimension_sizes, &
                      size(field_dimension_sizes), &
                      pelist(1), &
                      block = .false., &
                      tag=COMM_TAG_2)
        call mpp_sync_self(check=EVENT_RECV)
    endif

   !If the field_found flag is present, then return the value of the found
   !flag.  It is assumed that this value will be checked by the calling
   !routine.  If the field_found flag is not present and the field was not
   !found in the file, then throw a fatal error.
    if (present(field_found)) then
        field_found = found
    elseif (.not. found) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_get_field_size:" &
                       //" the inputted field "//trim(fieldname) &
                       //" was not found in the file "//trim(filename))
    endif

   !Deallocate local allocatables.
    deallocate(pelist)

    return
end subroutine fms_io_unstructured_get_field_size
