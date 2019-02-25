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
    MPP_TYPE_, dimension(:), intent(in) :: sbuf
    MPP_TYPE_, dimension(:), intent(inout) :: rbuf
    integer,   intent(in) :: scount, rcount

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call SYSTEM_CLOCK(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, MPP_TYPE_BYTELEN_)

end subroutine MPP_ALLTOALL_


subroutine MPP_ALLTOALLV_(sbuf, ssize, sdispl, rbuf, rsize, rdispl, pelist)
    MPP_TYPE_, intent(in) :: sbuf(:)
    MPP_TYPE_, intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call SYSTEM_CLOCK(start_tick)

    rbuf(:) = sbuf(:)

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

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call SYSTEM_CLOCK(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, MPP_TYPE_BYTELEN_)

end subroutine MPP_ALLTOALLW_
