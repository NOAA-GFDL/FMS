!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
#ifdef SYSTEM_CLOCK
#undef SYSTEM_CLOCK
#endif


program test   !test various aspects of mpp_mod
  use mpp_mod, only : mpp_init, mpp_exit, mpp_pe, mpp_npes, mpp_root_pe, mpp_commID, stdout
  use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_sync
  use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist, mpp_set_stack_size
  use mpp_mod, only : mpp_broadcast, mpp_transmit, mpp_sum, mpp_max, mpp_chksum, ALL_PES
  use mpp_mod, only : mpp_gather, mpp_error, FATAL, mpp_sync_self
  use platform_mod

  implicit none

  integer, parameter              :: n=1048576
  real, allocatable, dimension(:) :: a, b, c
  real, allocatable, dimension(:) :: d
  integer(i8_kind) :: locd
  integer                         :: tick, tick0, ticks_per_sec, id
  integer                         :: pe, npes, root, i, j, k, l, m, n2, istat
  integer                         :: out_unit
  real                            :: dt

  call mpp_init()
  call mpp_set_stack_size(3145746)
  pe = mpp_pe()
  npes = mpp_npes()
  root = mpp_root_pe()
  out_unit = stdout()

  call SYSTEM_CLOCK( count_rate=ticks_per_sec )

  if( pe.EQ.root ) print *, '------------------> Calling test_mpp_commID <------------------'
    call test_mpp_commID()
  if( pe.EQ.root ) print *, '------------------> Finished test_mpp_commID <------------------'

  if( pe.EQ.root ) print *, '------------------> Calling test_mpp_max <------------------'
    call test_mpp_max()
  if( pe.EQ.root ) print *, '------------------> Finished test_mpp_max <------------------'

!test of pointer sharing
  if( pe.EQ.root )then
      allocate( d(n) )
      locd = LOC(d)
  end if
  call mpp_broadcast(locd,root)
  if( pe.EQ.root )then
      call random_number(d)
  end if
  call mpp_sync()
!  call test_shared_pointers(locd,n)

  if( pe.EQ.root )then
      deallocate( d )
  end if

  call mpp_exit()

contains

  subroutine test_mpp_max

  allocate( a(n), b(n) )
  a = real(pe+1)
  print *, 'pe,     pe+1 =', pe, a(1)
  call mpp_max( a(1) )
  print *, 'pe, max(pe+1)=', pe, a(1)
  !pelist check
  call mpp_sync()
  flush(out_unit)
  if( npes.GE.2 )then
     if( pe.EQ.root )print *, 'Test of pelists: bcast, sum and max using PEs 0...npes-2 (excluding last PE)'
     call mpp_declare_pelist( (/(i,i=0,npes-2)/) )
     a = real(pe+1)
     if( pe.NE.npes-1 )call mpp_broadcast( a, n, npes-2, (/(i,i=0,npes-2)/) )
     print *, 'bcast(npes-1) from 0 to npes-2=', pe, a(1)
     a = real(pe+1)
     if( pe.NE.npes-1 )then
        call mpp_set_current_pelist( (/(i,i=0,npes-2)/) )
        id = mpp_clock_id( 'Partial mpp_sum' )
        call mpp_clock_begin(id)
        call mpp_sum( a(1:1000), 1000, (/(i,i=0,npes-2)/) )
        call mpp_clock_end  (id)
     end if
     if( pe.EQ.root )print *, 'sum(pe+1) from 0 to npes-2=', a(1)
     a = real(pe+1)
     if( pe.NE.npes-1 )call mpp_max( a(1), (/(i,i=0,npes-2)/) )
     if( pe.EQ.root )print *, 'max(pe+1) from 0 to npes-2=', a(1)
  end if
  call mpp_set_current_pelist()

  end subroutine test_mpp_max

  subroutine test_mpp_commID
    integer :: i, commID
    integer, allocatable :: pelist(:)

    ! Check if commID is valid (checks peset%id for default value, -1)
    if (mpp_commID().EQ.-1) then
       call mpp_error('test_mpp_commID', 'CommID Test after mpp_init Failed: mpp_commID returned an invalid ID', FATAL)
    end if

    ! Get commID for current pelist (should be global pelist) via and compare with mpp_commID
    allocate(pelist(npes))
    pelist = (/ (i, i=0, npes-1) /)
    call mpp_declare_pelist(pelist, commID = commID)

    if(mpp_commID().NE.commID) then
      call mpp_error('test_mpp_commID', 'Test failed: mpp_commID returned a different ID than expected.', FATAL)
    endif

    deallocate(pelist)

  end subroutine test_mpp_commID

  subroutine test_shared_pointers(locd,n)
    integer(i8_kind), intent(in) :: locd
    integer :: n
    real :: dd(n)
    pointer( p, dd )

    p = locd
    print *, 'TEST_SHARED_POINTERS: pe, locd=', pe, locd
!    print *, 'TEST_SHARED_POINTERS: pe, chksum(d)=', pe, mpp_chksum(dd,(/pe/))
    print *, 'TEST_SHARED_POINTERS: pe, sum(d)=', pe, sum(dd)
    return
  end subroutine test_shared_pointers
end program test
