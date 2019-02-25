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

subroutine MPP_ALLTOALL_(sbuf, scount, rbuf, rcount, pelist)
    MPP_TYPE_, intent(in) :: sbuf(:)
    MPP_TYPE_, intent(inout) :: rbuf(:)
    integer,   intent(in) :: scount, rcount
   
    integer, intent(in), optional :: pelist(0:)
    integer :: n

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    n = get_peset(pelist)

    if (current_clock .NE. 0) call SYSTEM_CLOCK(start_tick)

    if (verbose) call mpp_error(NOTE, 'MPP_ALLTOALL_: using MPI_Alltoall...')

    ! TODO: Message lengths greater than 1
    call MPI_Alltoall(sbuf, scount, MPI_TYPE_, rbuf, rcount, MPI_TYPE_, &
                      peset(n)%id, error)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, MPP_TYPE_BYTELEN_)

end subroutine MPP_ALLTOALL_


subroutine MPP_ALLTOALLV_(sbuf, ssize, sdispl, rbuf, rsize, rdispl, pelist)
    MPP_TYPE_, intent(in) :: sbuf(:)
    MPP_TYPE_, intent(inout) :: rbuf(:)

    ! TODO: Optionally set displacements to cumulative sums of ssize, rsize
    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)

    integer, intent(in), optional :: pelist(0:)
    integer :: n

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALLV_: You must first call mpp_init.')

    n = get_peset(pelist)

    if (current_clock .NE. 0) call SYSTEM_CLOCK(start_tick)

    if (verbose) call mpp_error(NOTE, 'MPP_ALLTOALLV_: using MPI_Alltoallv...')

    call MPI_Alltoallv(sbuf, ssize, sdispl, MPI_TYPE_, &
                       rbuf, rsize, rdispl, MPI_TYPE_, &
                       peset(n)%id, error)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, MPP_TYPE_BYTELEN_)

end subroutine MPP_ALLTOALLV_


subroutine MPP_ALLTOALLW_(sbuf, ssize, sdispl, stype, &
                          rbuf, rsize, rdispl, rtype, pelist)
    MPP_TYPE_, intent(in) :: sbuf(:)
    MPP_TYPE_, intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)
    type(mpp_type), intent(in) :: stype(:), rtype(:)
    integer, intent(in), optional :: pelist(0:)
    integer :: i, n

    integer, allocatable :: sendtypes(:), recvtypes(:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALLW_: You must first call mpp_init.')

    n = get_peset(pelist)

    if (current_clock .NE. 0) call SYSTEM_CLOCK(start_tick)

    if (verbose) call mpp_error(NOTE, 'MPP_ALLTOALLW_: using MPI_Alltoallw...')

    ! Convert mpp_types to MPI datatype IDs
    ! NOTE: sendtypes and recvtypes must be the same size
    allocate(sendtypes(size(stype)))
    allocate(recvtypes(size(rtype)))
    do i = 1, size(stype)
        sendtypes(i) = stype(i)%id
        recvtypes(i) = rtype(i)%id
    end do

    call MPI_Alltoallw(sbuf, ssize, sdispl, sendtypes, &
                       rbuf, rsize, rdispl, recvtypes, &
                       peset(n)%id, error)

    deallocate(sendtypes, recvtypes)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, MPP_TYPE_BYTELEN_)

end subroutine MPP_ALLTOALLW_
