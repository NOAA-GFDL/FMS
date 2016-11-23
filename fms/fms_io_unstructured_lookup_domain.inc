!----------
!ug support

!>This routine returns the index of the inputted unstructured domain in
!!the domain_UG_array module array, or else the value -1.
function fms_io_unstructured_lookup_domain(domain) &
                                           result(domain_index)

   !Inputs/Outputs
    type(domainUG),intent(in) :: domain       !<An unstructured mpp domain.
    integer(INT_KIND)         :: domain_index !<The index of the domain in the domain_UG_array module array.

   !Local variables
    integer(INT_KIND) :: i !<Loop variable.

   !domain_UG_array_size is a module variable.
   !domain_UG_array is a module variable.
    domain_index = -1
    do i = 1,domain_UG_array_size
        if (mpp_compare_UG_domains(domain,domain_UG_array(i))) then
            domain_index = i
            exit
        endif
    enddo

    return
end function fms_io_unstructured_lookup_domain

!----------
