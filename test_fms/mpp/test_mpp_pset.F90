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


program test_mpp_pset
  use mpp_mod, only: mpp_init, mpp_exit, mpp_pe, mpp_npes, stderr, stdout, &
       mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_error, FATAL
 use mpp_mod, only : input_nml_file
  use mpp_pset_mod !, only: mpp_pset_type, mpp_pset_create, mpp_pset_root, &
!       mpp_pset_broadcast_ptr, mpp_pset_segment_array, mpp_pset_sync, &
!       mpp_pset_stack_push, mpp_pset_print_chksum, mpp_pset_delete
  implicit none
!test program demonstrates how to create PSETs
!  how to distribute allocatable arrays
!  how to distribute automatic arrays
  integer, parameter :: n=96 !divisible by lots of numbers
  real, allocatable, dimension(:,:,:) :: a, b, cc
  real :: c(n,n,n)
  logical :: opened
#ifdef use_CRI_pointers
  pointer( ptr_c, c )
#endif
  integer, pointer :: ptr !useless declaration, but it will compile
  integer :: i, j, k, ks, ke
!MPP
  integer :: pe, npes
!MPP_PSET
  type(mpp_pset_type) :: pset
  logical :: root
!clocks
  integer :: id_full, id_alloc, id_auto, id
  integer :: out_unit, errunit, io_status
  integer :: test_number
  integer :: unit=7

namelist / test_mpp_pset_nml / test_number

  call mpp_init()

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, test_mpp_pset_nml, iostat=io_status)
#else
  do
     inquire( unit=unit, opened=opened )
     if( .NOT.opened )exit
     unit = unit + 1
     if( unit.EQ.100 )call mpp_error( FATAL, 'Unable to locate unit number.' )
  end do
  open( unit=unit, file='input.nml', iostat=io_status )
  read( unit,test_mpp_pset_nml, iostat=io_status )
  close(unit)
#endif

      if (io_status > 0) then
         call mpp_error(FATAL,'=>test_mpp_domains: Error reading input.nml')
endif

  pe = mpp_pe()
  npes = mpp_npes()
  out_unit = stdout()
  errunit = stderr()
  write( out_unit,'(a,i6)' )'Starting MPP_PSET unit test, npes=', npes
  call mpp_pset_create( npes, pset )
  root = mpp_pset_root(pset)
  id_full = mpp_clock_id( 'Full array' )
  id_alloc = mpp_clock_id( 'Allocatable array, PSETs' )
  id_auto = mpp_clock_id( 'Automatic array, PSETs' )
!allocate a and b
  allocate( a(n,n,n) )
  allocate( b(n,n,n) )

!allocate shared array c
  if( root )then
      allocate( cc(n,n,n) )
#ifdef use_CRI_pointers
      ptr = LOC(cc)
#endif
end if

!  call mpp_pset_broadcast_ptr( pset, ptr )

#ifdef use_CRI_pointers
  ptr_c = ptr
#endif

!initialize a and b
  call RANDOM_NUMBER(a)
  call mpp_clock_begin(id_full)
  do k = 1,n
     do j = 1,n
        do i = 1,n
           b(i,j,k) = 2*a(i,j,k)
        end do
     end do
  end do
  call mpp_clock_end(id_full)

   if (test_number == 1) then
      ! Testing how to distribute allocatable arrays
      id = id_alloc
   else if (test_number == 2) then
          ! Testing how you create shared auto arrays
#ifdef use_CRI_pointers
    pointer( pd, c )
    call mpp_pset_stack_push( pset, pd, size(c) )
#endif
          id = id_auto
   endif

!divide up among PSETs
  call mpp_pset_segment_array( pset, 1, n, ks, ke )
  write( errunit,'(a,4i6)' )'pe, n, ks, ke=', pe, n, ks, ke
  call mpp_clock_begin(id)
  do k = ks,ke
     do j = 1,n
        do i = 1,n
           c(i,j,k) = 2*a(i,j,k)
        end do
     end do
  end do
  call mpp_pset_sync(pset)
  call mpp_clock_end(id)

  write( errunit,'(a,i6,2es23.15)' )'b, c should be equal: pe b c=', &
       pe, sum(b), sum(c)
  call mpp_pset_print_chksum( pset, 'test_alloc', c(:,:,ks:ke) )
  call mpp_pset_delete(pset)
  call mpp_exit()

end program test_mpp_pset
