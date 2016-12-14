!----------
!ug support

!>Find the file unit for an inputted file, searching for its variants.  If the
!!file is not found, then throw a fatal error.
subroutine fms_io_unstructured_file_unit(filename, &
                                         funit, &
                                         domain)

   !Inputs/Outputs
    character(len=*),intent(in)   :: filename   !<The name of a file.
    integer(INT_KIND),intent(out) :: funit      !<The file unit for the input file.
    type(domainUG),intent(in)     :: domain     !<An unstructured mpp domain.

   !Local variables
    logical(INT_KIND)  :: found_file  !<Flag telling if the file exists.
    character(len=256) :: actual_file !<Name of the found file.
    logical(INT_KIND)  :: read_dist   !<Flag telling if the file is "distributed" (has IO domain tile id appended to the end).
    integer(INT_KIND)  :: nfile       !<Index of the inputted file in the "files_read" module array.

   !Get the actual name of the file, searching for all possible variants. If
   !the file is not found, then throw a fatal error.
    found_file = fms_io_unstructured_get_file_name(filename, &
                                                   actual_file, &
                                                   read_dist, &
                                                   domain)
    if (.not. found_file) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_file_unit: file " &
                       //trim(filename)//" and corresponding variants are" &
                       //" not found.")
    endif

   !Get the file unit for the input file.
    call fms_io_unstructured_get_file_unit(actual_file, &
                                           funit, &
                                           nfile, &
                                           read_dist, &
                                           domain)

    return
end subroutine fms_io_unstructured_file_unit

!----------
