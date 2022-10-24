! -*-f90-*-

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

!> @file
!> @brief MPI implementation of @ref mpp_alltoall routines

!> @addtogroup mpp_mod
!> @{

!> Wrapper for mpi_alltoall routine, sends data from all to all processes
subroutine MPP_ALLTOALL_(sbuf, scount, rbuf, rcount, pelist)
    MPP_TYPE_, intent(in) :: sbuf(:) !< data to send
    MPP_TYPE_, intent(inout) :: rbuf(:) !< received data
    integer,   intent(in) :: scount !< length of buffer data to send from each process
    integer,   intent(in) :: rcount !< length of buffer data to recieve from each proces

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

!> Wrapper for mpi_alltoallv, sends data from all to all processes with vector displacement
subroutine MPP_ALLTOALLV_(sbuf, ssize, sdispl, rbuf, rsize, rdispl, pelist)
    MPP_TYPE_, intent(in) :: sbuf(:) !< data to send
    MPP_TYPE_, intent(inout) :: rbuf(:) !< received data

    ! TODO: Optionally set displacements to cumulative sums of ssize, rsize
    integer, intent(in) :: ssize(:) !< array containing number of elements to send to each respective process
    integer, intent(in) :: rsize(:) !< array containing number of elements to recieve from each respective process
    integer, intent(in) :: sdispl(:)!< displacement of data sent to each respective process
    integer, intent(in) :: rdispl(:)!< displacement of data received from each respective process

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

!> Wrapper for mpi_alltoallw, sends data from all to all processes with given data types,
!! displacements and block sizes
subroutine MPP_ALLTOALLW_(sbuf, ssize, sdispl, stype, &
                          rbuf, rsize, rdispl, rtype, pelist)
    MPP_TYPE_, intent(in) :: sbuf(:) !< data to send
    MPP_TYPE_, intent(inout) :: rbuf(:) !< recieved data

    integer, intent(in) :: ssize(:) !< array containing number of elements to send to each respective process
    integer, intent(in) :: rsize(:) !< array containing number of elements to recieve from each respective process
    integer, intent(in) :: sdispl(:)!< displacement of data sent to each respective process
    integer, intent(in) :: rdispl(:)!< displacement of data received from each respective process
    type(mpp_type), intent(in) :: stype(:) !< mpp data types to send to each respective process
    type(mpp_type), intent(in) :: rtype(:) !< mpp data types to recieve from each respective process
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
!> @}
