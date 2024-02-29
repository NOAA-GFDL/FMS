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
!> @author Eric Stofferahn
!> @brief Test mpp_min and mpp_max functions for various precisions of
!! reals
program test

  use mpp_mod, only : mpp_init, mpp_pe, mpp_npes, mpp_root_pe, stdout
  use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_sync
  use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist, mpp_set_stack_size
  use mpp_mod, only : mpp_broadcast, mpp_sum, mpp_min, mpp_max
  use mpp_mod, only : mpp_error, FATAL
  use platform_mod

  implicit none

  integer, parameter              :: n=1048576
  real(kind=r4_kind), allocatable, dimension(:) :: a4
  real(kind=r8_kind), allocatable, dimension(:) :: a8
  integer(kind=i4_kind), allocatable, dimension(:) :: b4
  integer(kind=i8_kind), allocatable, dimension(:) :: b8
  integer                         :: id, pe, npes, root, i, out_unit, ierr

  call mpp_init(0)
  call mpp_set_stack_size(3145746)
  pe = mpp_pe()
  npes = mpp_npes()
  root = mpp_root_pe()
  out_unit = stdout()
  allocate( a4(n), a8(n), b4(n), b8(n) )

  if( pe.EQ.root ) print *, '-> Calling test_mpp_max_r4 <-------------------'
    call test_mpp_max_r4()
  if( pe.EQ.root ) print *, '-> test_mpp_max_r4: <------------------ Passed!'

  if( npes.GE.2 ) then
    if( pe.EQ.root ) print *, '-> Calling test_mpp_max_with_pe_r4 <-----------'
      call test_mpp_max_with_pe_r4()
    if( pe.EQ.root ) print *, '-> test_mpp_max_with_pe_r4: <---------- Passed!'
  else
    if( pe.EQ.root ) print *, '-> test_mpp_max_with_pe_r4: <- (one pe) Skipped'
  end if

  if( pe.EQ.root ) print *, '-> Calling test_mpp_max_r8 <-------------------'
    call test_mpp_max_r8()
  if( pe.EQ.root ) print *, '-> test_mpp_max_r8: <------------------ Passed!'

  if( npes.GE.2 ) then
    if( pe.EQ.root ) print *, '-> Calling test_mpp_max_with_pe_r8 <-----------'
      call test_mpp_max_with_pe_r8()
    if( pe.EQ.root ) print *, '-> test_mpp_max_with_pe_r8: <---------- Passed!'
  else
    if( pe.EQ.root ) print *, '-> test_mpp_max_with_pe_r8: <- (one pe) Skipped'
  end if

  if( pe.EQ.root ) print *, '-> Calling test_mpp_min_r4 <-------------------'
    call test_mpp_min_r4()
  if( pe.EQ.root ) print *, '-> test_mpp_min_r4: <------------------ Passed!'

  if( npes.GE.2 ) then
    if( pe.EQ.root ) print *, '-> Calling test_mpp_min_with_pe_r4 <-----------'
      call test_mpp_min_with_pe_r4()
    if( pe.EQ.root ) print *, '-> test_mpp_min_with_pe_r4: <---------- Passed!'
  else
    if( pe.EQ.root ) print *, '-> test_mpp_min_with_pe_r4: <- (one pe) Skipped'
  end if

  if( pe.EQ.root ) print *, '-> Calling test_mpp_min_r8 <-------------------'
    call test_mpp_min_r8()
  if( pe.EQ.root ) print *, '-> test_mpp_min_r8: <------------------ Passed!'

  if( npes.GE.2 ) then
    if( pe.EQ.root ) print *, '-> Calling test_mpp_min_with_pe_r8 <-----------'
      call test_mpp_min_with_pe_r8()
    if( pe.EQ.root ) print *, '-> test_mpp_min_with_pe_r8: <---------- Passed!'
  else
    if( pe.EQ.root ) print *, '-> test_mpp_min_with_pe_r8: <- (one pe) Skipped'
  end if

  if( pe.EQ.root ) print *, '-> Calling test_mpp_max_i4 <-------------------'
    call test_mpp_max_i4()
  if( pe.EQ.root ) print *, '-> test_mpp_max_i4: <------------------ Passed!'

  if( npes.GE.2 ) then
    if( pe.EQ.root ) print *, '-> Calling test_mpp_max_with_pe_i4 <-----------'
      call test_mpp_max_with_pe_i4()
    if( pe.EQ.root ) print *, '-> test_mpp_max_with_pe_i4: <---------- Passed!'
  else
    if( pe.EQ.root ) print *, '-> test_mpp_max_with_pe_i4: <- (one pe) Skipped'
  end if

  if( pe.EQ.root ) print *, '-> Calling test_mpp_max_i8 <-------------------'
    call test_mpp_max_i8()
  if( pe.EQ.root ) print *, '-> test_mpp_max_i8: <------------------ Passed!'

  if( npes.GE.2 ) then
    if( pe.EQ.root ) print *, '-> Calling test_mpp_max_with_pe_i8 <-----------'
      call test_mpp_max_with_pe_i8()
    if( pe.EQ.root ) print *, '-> test_mpp_max_with_pe_i8: <---------- Passed!'
  else
    if( pe.EQ.root ) print *, '-> test_mpp_max_with_pe_i8: <- (one pe) Skipped'
  end if

  if( pe.EQ.root ) print *, '-> Calling test_mpp_min_i4 <-------------------'
    call test_mpp_min_i4()
  if( pe.EQ.root ) print *, '-> test_mpp_min_i4: <------------------ Passed!'

  if( npes.GE.2 ) then
    if( pe.EQ.root ) print *, '-> Calling test_mpp_min_with_pe_i4 <-----------'
      call test_mpp_min_with_pe_i4()
    if( pe.EQ.root ) print *, '-> test_mpp_min_with_pe_i4: <---------- Passed!'
  else
    if( pe.EQ.root ) print *, '-> test_mpp_min_with_pe_i4: <- (one pe) Skipped'
  end if

  if( pe.EQ.root ) print *, '-> Calling test_mpp_min_i8 <-------------------'
    call test_mpp_min_i8()
  if( pe.EQ.root ) print *, '-> test_mpp_min_i8: <------------------ Passed!'

  if( npes.GE.2 ) then
    if( pe.EQ.root ) print *, '-> Calling test_mpp_min_with_pe_i8 <-----------'
      call test_mpp_min_with_pe_i8()
    if( pe.EQ.root ) print *, '-> test_mpp_min_with_pe_i8: <---------- Passed!'
  else
    if( pe.EQ.root ) print *, '-> test_mpp_min_with_pe_i8: <- (one pe) Skipped'
  end if

  deallocate( a4, a8, b4, b8 )
  call MPI_FINALIZE(ierr)

contains

  subroutine test_mpp_max_r4
    a4 = real(pe+1, kind=r4_kind)
    call mpp_max( a4(1) )
    if (a4(1).NE.real(npes, kind=r4_kind)) then
      call mpp_error(FATAL, "The r4 mpp_max function for all npes did not return the appropriate answer")
    end if
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_max_r4

  subroutine test_mpp_max_with_pe_r4
    call mpp_declare_pelist( (/(i,i=0,npes-2)/) )
    if(pe.NE.npes-1) call mpp_set_current_pelist( (/(i,i=0,npes-2)/) )
    a4 = real(pe+1, kind=r4_kind)
    if( pe.NE.npes-1 ) then
      call mpp_max( a4(1), (/(i,i=0,npes-2)/) )
      if (a4(1).NE.real(npes-1, kind=r4_kind)) then
        call mpp_error(FATAL, "The r4 mpp_max function for all but the last pe did not return the appropriate answer")
      end if
    end if
    call mpp_set_current_pelist()
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_max_with_pe_r4

  subroutine test_mpp_max_r8
    a8 = real(pe+1, kind=r8_kind)
    call mpp_max( a8(1) )
    if (a8(1).NE.real(npes, kind=r8_kind)) then
      call mpp_error(FATAL, "The r8 mpp_max function for all npes did not return the appropriate answer")
    end if
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_max_r8

  subroutine test_mpp_max_with_pe_r8
    call mpp_declare_pelist( (/(i,i=0,npes-2)/) )
    if(pe.NE.npes-1) call mpp_set_current_pelist( (/(i,i=0,npes-2)/) )
    a8 = real(pe+1, kind=r8_kind)
    if( pe.NE.npes-1 ) then
      call mpp_max( a8(1), (/(i,i=0,npes-2)/) )
      if (a8(1).NE.real(npes-1, kind=r8_kind)) then
        call mpp_error(FATAL, "The r8 mpp_max function for all but the last pe did not return the appropriate answer")
      end if
    end if
    call mpp_set_current_pelist()
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_max_with_pe_r8

  subroutine test_mpp_min_r4
    a4 = real(pe+1, kind=r4_kind)
    call mpp_min( a4(1) )
    if (a4(1).NE.real(1, kind=r4_kind)) then
      call mpp_error(FATAL, "The r4 mpp_min function for all npes did not return the appropriate answer")
    end if
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_min_r4

  subroutine test_mpp_min_with_pe_r4
    call mpp_declare_pelist( (/(i,i=0,npes-2)/) )
    if(pe.NE.npes-1) call mpp_set_current_pelist( (/(i,i=0,npes-2)/) )
    a4 = real(pe+1, kind=r4_kind)
    if( pe.NE.npes-1 ) then
      call mpp_min( a4(1), (/(i,i=0,npes-2)/) )
      if (a4(1).NE.real(1, kind=r4_kind)) then
        call mpp_error(FATAL, "The r4 mpp_min function for all but the last pe did not return the appropriate answer")
      end if
    end if
    call mpp_set_current_pelist()
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_min_with_pe_r4

  subroutine test_mpp_min_r8
    a8 = real(pe+1, kind=r8_kind)
    call mpp_min( a8(1) )
    if (a8(1).NE.real(1, kind=r8_kind)) then
      call mpp_error(FATAL, "The r8 mpp_min function for all npes did not return the appropriate answer")
    end if
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_min_r8

  subroutine test_mpp_min_with_pe_r8
    call mpp_declare_pelist( (/(i,i=0,npes-2)/) )
    if(pe.NE.npes-1) call mpp_set_current_pelist( (/(i,i=0,npes-2)/) )
    a8 = real(pe+1, kind=r8_kind)
    if( pe.NE.npes-1 ) then
      call mpp_min( a8(1), (/(i,i=0,npes-2)/) )
      if (a8(1).NE.real(1, kind=r8_kind)) then
        call mpp_error(FATAL, "The r8 mpp_min function for all but the last pe did not return the appropriate answer")
      end if
    end if
    call mpp_set_current_pelist()
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_min_with_pe_r8

  subroutine test_mpp_max_i4
    b4 = int(pe+1, kind=i4_kind)
    call mpp_max( b4(1) )
    if (b4(1).NE.int(npes, kind=i4_kind)) then
      call mpp_error(FATAL, "The i4 mpp_max function for all npes did not return the appropriate answer")
    end if
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_max_i4

  subroutine test_mpp_max_with_pe_i4
    call mpp_declare_pelist( (/(i,i=0,npes-2)/) )
    if(pe.NE.npes-1) call mpp_set_current_pelist( (/(i,i=0,npes-2)/) )
    b4 = int(pe+1, kind=i4_kind)
    if( pe.NE.npes-1 ) then
      call mpp_max( b4(1), (/(i,i=0,npes-2)/) )
      if (b4(1).NE.int(npes-1, kind=i4_kind)) then
        call mpp_error(FATAL, "The i4 mpp_max function for all but the last pe did not return the appropriate answer")
      end if
    end if
    call mpp_set_current_pelist()
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_max_with_pe_i4

  subroutine test_mpp_max_i8
    b8 = int(pe+1, kind=i8_kind)
    call mpp_max( b8(1) )
    if (b8(1).NE.int(npes, kind=i8_kind)) then
      call mpp_error(FATAL, "The i8 mpp_max function for all npes did not return the appropriate answer")
    end if
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_max_i8

  subroutine test_mpp_max_with_pe_i8
    call mpp_declare_pelist( (/(i,i=0,npes-2)/) )
    if(pe.NE.npes-1) call mpp_set_current_pelist( (/(i,i=0,npes-2)/) )
    b8 = int(pe+1, kind=i8_kind)
    if( pe.NE.npes-1 ) then
      call mpp_max( b8(1), (/(i,i=0,npes-2)/) )
      if (b8(1).NE.int(npes-1, kind=i8_kind)) then
        call mpp_error(FATAL, "The i8 mpp_max function for all but the last pe did not return the appropriate answer")
      end if
    end if
    call mpp_set_current_pelist()
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_max_with_pe_i8

  subroutine test_mpp_min_i4
    b4 = int(pe+1, kind=i4_kind)
    call mpp_min( b4(1) )
    if (b4(1).NE.int(1, kind=i4_kind)) then
      call mpp_error(FATAL, "The i4 mpp_min function for all npes did not return the appropriate answer")
    end if
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_min_i4

  subroutine test_mpp_min_with_pe_i4
    call mpp_declare_pelist( (/(i,i=0,npes-2)/) )
    if(pe.NE.npes-1) call mpp_set_current_pelist( (/(i,i=0,npes-2)/) )
    b4 = int(pe+1, kind=i4_kind)
    if( pe.NE.npes-1 ) then
      call mpp_min( b4(1), (/(i,i=0,npes-2)/) )
      if (b4(1).NE.int(1, kind=i4_kind)) then
        call mpp_error(FATAL, "The i4 mpp_min function for all but the last pe did not return the appropriate answer")
      end if
    end if
    call mpp_set_current_pelist()
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_min_with_pe_i4

  subroutine test_mpp_min_i8
    b8 = int(pe+1, kind=i8_kind)
    call mpp_min( b8(1) )
    if (b8(1).NE.int(1, kind=i8_kind)) then
      call mpp_error(FATAL, "The i8 mpp_min function for all npes did not return the appropriate answer")
    end if
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_min_i8

  subroutine test_mpp_min_with_pe_i8
    call mpp_declare_pelist( (/(i,i=0,npes-2)/) )
    if(pe.NE.npes-1) call mpp_set_current_pelist( (/(i,i=0,npes-2)/) )
    b8 = int(pe+1, kind=i8_kind)
    if( pe.NE.npes-1 ) then
      call mpp_min( b8(1), (/(i,i=0,npes-2)/) )
      if (b8(1).NE.int(1, kind=i8_kind)) then
        call mpp_error(FATAL, "The i8 mpp_min function for all but the last pe did not return the appropriate answer")
      end if
    end if
    call mpp_set_current_pelist()
    call mpp_sync()
    call flush(out_unit)
  end subroutine test_mpp_min_with_pe_i8

end program test
