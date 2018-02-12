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

!>Return the file unit and index in the "files_read" module array for the
!!inputted file.  If the file does not currently exist in the "files_read"
!!array (i.e., it is not currenly open), then open it.
subroutine fms_io_unstructured_get_file_unit(filename, &
                                             funit, &
                                             index_file, &
                                             read_dist, &
                                             domain)

   !Inputs/Outputs
    character(len=*),intent(in)   :: filename   !<Name of the file to be read from.
    integer(INT_KIND),intent(out) :: funit      !<File unit for the inputted file.
    integer(INT_KIND),intent(out) :: index_file !<Index of the inputted file in the "files_read" module array.
    logical(INT_KIND),intent(in)  :: read_dist  !<Flag telling if the IO domain tile id string exists at the end of the inputted file name.
    type(domainUG),intent(in)     :: domain     !<An unstructured mpp domain.

   !Local variables
    integer(INT_KIND) :: i !<Loop variable.

   !Check if the file exists in the "files_read" module array.  If the file
   !is found in the array, this implies that the file was opened at some
   !point.  If the file is still open, then return.  If not, then throw a
   !fatal error.
   !num_files_r is a module variable.
   !files_read is a module variable.
    do i = 1,num_files_r
        if (trim(files_read(i)%name) .eq. trim(filename)) then
            index_file = i
            funit = files_read(index_file)%unit
            if (.not. mpp_file_is_opened(funit)) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_get_file_unit:" &
                               //" the file "//trim(filename) &
                               //" was previously opened, but is not" &
                               //" currently open.")
            endif
            return
        endif
    enddo

   !Make sure that the number of files open for reading will not exceed the
   !maximum number allowed.
   !max_files_r is a module variable.
    if (num_files_r .eq. max_files_r) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_get_file_unit: the number of" &
                       //" files currently open for reading exceeds" &
                       //" max_files_r.  Pleaes increase this value via" &
                       //" the fms_io_nml namelist.")
    endif

   !Store the filename and initialize an array for the file fields.
    num_files_r = num_files_r + 1
    files_read(num_files_r)%name = trim(filename)
    allocate(files_read(num_files_r)%var(max_fields))
    files_read(num_files_r)%nvar = 0
    index_file = num_files_r

   !Open the file.
    if (read_dist) then
        call mpp_open(funit, &
                      trim(filename), &
                      form=form, &
                      action=MPP_RDONLY, &
                      threading=MPP_MULTI, &
                      fileset=MPP_MULTI, &
                      domain_ug=domain)
    else
        call mpp_open(funit, &
                      trim(filename), &
                      form=form, &
                      action=MPP_RDONLY, &
                      threading=MPP_MULTI, &
                      fileset=MPP_SINGLE)
    endif

   !Store the file unit returned by mpp_open.
    files_read(index_file)%unit = funit

    return
end subroutine fms_io_unstructured_get_file_unit

!----------
