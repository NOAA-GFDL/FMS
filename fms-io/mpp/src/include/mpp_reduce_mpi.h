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
    subroutine MPP_REDUCE_0D_( a, pelist )
!find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast to all PEs
      MPP_TYPE_, intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      MPP_TYPE_ :: work

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_0D: You must first call mpp_init.' )
      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( debug .and. (current_clock.NE.0) )call SYSTEM_CLOCK(start_tick)
      if( verbose )call mpp_error( NOTE, 'MPP_REDUCE_0D_: using MPI_ALLREDUCE...' )
      call MPI_ALLREDUCE( a, work, 1, MPI_TYPE_, MPI_REDUCE_, peset(n)%id, error )
      a = work
      if( debug .and. (current_clock.NE.0) )call increment_current_clock( EVENT_ALLREDUCE, MPP_TYPE_BYTELEN_ )
      return
    end subroutine MPP_REDUCE_0D_

    subroutine MPP_REDUCE_1D_( a, length, pelist )
!find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast to all PEs
      MPP_TYPE_, intent(inout) :: a(:)
      integer,   intent(in)    :: length
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      MPP_TYPE_ :: work(length)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_1D: You must first call mpp_init.' )
      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( debug .and. (current_clock.NE.0) )call SYSTEM_CLOCK(start_tick)
      if( verbose )call mpp_error( NOTE, 'MPP_REDUCE_1D_: using MPI_ALLREDUCE...' )
      call MPI_ALLREDUCE( a, work, length, MPI_TYPE_, MPI_REDUCE_, peset(n)%id, error )
      a(1:length) = work(1:length)
      if( debug .and. (current_clock.NE.0) )call increment_current_clock( EVENT_ALLREDUCE, MPP_TYPE_BYTELEN_ )
      return
    end subroutine MPP_REDUCE_1D_


