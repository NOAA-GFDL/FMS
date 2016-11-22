!----------
!ug support

!>For a file associated with an unstructured mpp domain, check if the file
!!or any of its variants exist.
function fms_io_unstructured_file_exist(file_name, &
                                        domain) &
                                        result(does_file_exist)

   !Inputs/Outputs
    character(len=*),intent(in) :: file_name       !<The name of a file.
    type(domainUG),intent(in)   :: domain          !<An unstructured mpp domain.
    logical(INT_KIND)           :: does_file_exist !<Flag telling if the inputted file or any of its variants exist.

   !Local variables
    character(len=256) :: fname     !<Actual name of the found file.
    logical(INT_KIND)  :: read_dist !<Flag that tells if the file is "distributed" (has I/O domain tile id appended onto the end).

   !Check if the file exists.
    does_file_exist = fms_io_unstructured_get_file_name(file_name, &
                                                        fname, &
                                                        read_dist, &
                                                        domain)

    return
end function fms_io_unstructured_file_exist

!----------
