!----------
!ug support

!------------------------------------------------------------------------------
function mpp_domain_UG_is_tile_root_pe(domain) result(is_root)

   !<Inputs/Outputs
    type(domainUG),intent(in) :: domain
    logical(INT_KIND)         :: is_root

    if (domain%pe .eq. domain%tile_root_pe) then
        is_root = .true.
    else
        is_root = .false.
    endif

    return
end function mpp_domain_UG_is_tile_root_pe

!------------------------------------------------------------------------------
!HELP: There needs to be a subroutine to return the "io_layout" for
!      an unstructured domain, so I made one.  Someone should check
!      to see if this is correct.
function mpp_get_io_domain_UG_layout(domain) result(io_layout)

   !<Inputs/Outputs
    type(domainUG),intent(in) :: domain
    integer(INT_KIND)         :: io_layout

    io_layout = domain%io_layout

    return
end function

!------------------------------------------------------------------------------
!HELP: There needs to be a subroutine to deallocate an unstructured grid,
!       so I made one.  Someone should check to see if this is correct.
subroutine mpp_deallocate_domainUG(domain)

   !<Inputs/Outputs
    type(domainUG),intent(inout) :: domain

   !<Local variables
    integer(INT_KIND) :: i !<Loop variable.

    if (associated(domain%list)) then
        deallocate(domain%list)
        domain%list => null()
    endif

    if (associated(domain%io_domain)) then
        if (associated(domain%io_domain%list)) then
            deallocate(domain%io_domain%list)
            domain%io_domain%list => null()
        endif
        deallocate(domain%io_domain)
        domain%io_domain => null()
    endif

    if (associated(domain%UG2SG%recv)) then
        domain%UG2SG%recv => null()
    endif
    if (associated(domain%UG2SG%send)) then
        domain%UG2SG%send => null()
    endif

    if (associated(domain%SG2UG%recv)) then
        do i = 1,size(domain%SG2UG%recv)
            if (associated(domain%SG2UG%recv(i)%i)) then
                deallocate(domain%SG2UG%recv(i)%i)
                domain%SG2UG%recv(i)%i => null()
            endif
            if (associated(domain%SG2UG%recv(i)%j)) then
                deallocate(domain%SG2UG%recv(i)%j)
                domain%SG2UG%recv(i)%j => null()
            endif
        enddo
        deallocate(domain%SG2UG%recv)
        domain%SG2UG%recv => null()
    endif
    if (associated(domain%SG2UG%send)) then
        do i = 1,size(domain%SG2UG%send)
            if (associated(domain%SG2UG%send(i)%i)) then
                deallocate(domain%SG2UG%send(i)%i)
                domain%SG2UG%send(i)%i => null()
            endif
            if (associated(domain%SG2UG%send(i)%j)) then
                deallocate(domain%SG2UG%send(i)%j)
                domain%SG2UG%send(i)%j => null()
            endif
        enddo
        deallocate(domain%SG2UG%send)
        domain%SG2UG%send => null()
    endif

    if (associated(domain%grid_index)) then
        deallocate(domain%grid_index)
        domain%grid_index => null()
    endif

    if (associated(domain%SG_domain)) then
        domain%SG_domain => null()
    endif

    return
end subroutine mpp_deallocate_domainUG

!------------------------------------------------------------------------------
!>HELP: I need a function that checks if 2 different unstructured domains are
!!      equal, so I made one.  Someone should check this.
function mpp_compare_UG_domains(domain1, &
                                domain2) &
                                result(domains_are_equal)

   !Inputs/Outputs
    type(domainUG),intent(in) :: domain1 !<An unstructured mpp domain.
    type(domainUG),intent(in) :: domain2 !<An unstructured mpp domain.
    logical(INT_KIND)         :: domains_are_equal

    if (associated(domain1%SG_domain) .and. associated(domain2%SG_domain)) then
        if (domain1%SG_domain .ne. domain2%SG_domain) then
            domains_are_equal = .false.
            return
        endif
    elseif (associated(domain1%SG_domain) .and. .not. &
            associated(domain2%SG_domain)) then
        domains_are_equal = .false.
        return
    elseif (.not. associated(domain1%SG_domain) .and. &
            associated(domain2%SG_domain)) then
        domains_are_equal = .false.
        return
    endif

    if (associated(domain1%list) .and. associated(domain2%list)) then
        if (size(domain1%list) .ne. size(domain2%list)) then
            domains_are_equal = .false.
            return
        endif
    elseif (associated(domain1%list) .and. .not. associated(domain2%list)) then
        domains_are_equal = .false.
        return
    elseif (.not. associated(domain1%list) .and. associated(domain2%list)) then
        domains_are_equal = .false.
        return
    endif

    if (associated(domain1%grid_index) .and. associated(domain2%grid_index)) then
        if (size(domain1%grid_index) .ne. size(domain2%grid_index)) then
            domains_are_equal = .false.
            return
        endif
    elseif (associated(domain1%grid_index) .and. .not. associated(domain2%grid_index)) then
        domains_are_equal = .false.
        return
    elseif (.not. associated(domain1%grid_index) .and. associated(domain2%grid_index)) then
        domains_are_equal = .false.
        return
    endif

    if (domain1%pe .ne. domain2%pe .or. &
        domain1%pos .ne. domain2%pos .or. &
        domain1%ntiles .ne. domain2%ntiles .or. &
        domain1%tile_id .ne. domain2%tile_id .or. &
        domain1%tile_root_pe .ne. domain2%tile_root_pe) then
        domains_are_equal = .false.
        return
    endif

    domains_are_equal = .true.

    return
end function mpp_compare_UG_domains

!------------------------------------------------------------------------------

!----------
