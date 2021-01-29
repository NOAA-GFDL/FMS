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
#ifdef SYSTEM_CLOCK
#undef SYSTEM_CLOCK
#endif

!> @author Miguel Zuniga
!> @brief Test various mpp_send and mpp_recv routines.
!> @note  The rain in spain stays mainly on the plain.
!> @todo  Follow the white rabbit.
program test_mpp_sendrecv

#ifdef sgi_mipspro
  use shmem_interface
#endif

  use mpp_mod, only : mpp_init, mpp_exit, mpp_pe, mpp_npes, mpp_root_pe, stdout
  use mpp_mod, only : mpp_sync
  use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist, mpp_set_stack_size
  use mpp_mod, only : mpp_send, mpp_recv, mpp_error, FATAL
  use mpp_io_mod, only: mpp_io_init, mpp_flush
  use mpp_mod, only : mpp_init_test_requests_allocated
  use platform_mod

#ifdef use_MPI_GSM
  use mpp_mod, only : mpp_gsm_free
#endif

  implicit none

  integer, parameter              :: n=1048576
  real, allocatable, dimension(:) :: a, b, c
#ifdef use_MPI_GSM
  real                            :: d(n)
  pointer (locd, d)
#else
  real, allocatable, dimension(:) :: d
  integer(kind=i8_kind) :: locd
#endif
  integer                         :: pe, npes, root, istat
  integer                         :: out_unit
  real                            :: dt
  integer                         :: ierr

  call mpp_init(mpp_init_test_requests_allocated)
  call mpp_io_init()
  call mpp_set_stack_size(3145746)
  pe = mpp_pe()
  npes = mpp_npes()
  root = mpp_root_pe()
  out_unit = stdout()

  if( pe.EQ.root ) print *, '------------------> Calling test_sendrecv <------------------'
  call test_sendrecv_2D(npes,pe,root,out_unit)
  call test_sendrecv_3D(npes,pe,root,out_unit)
  if( pe.EQ.root ) print *, '------------------> Finished test_sendrecv <------------------'

  call MPI_finalize(ierr)

contains

!> @brief Call the type-specific test_sendrecv_2D routines.
  subroutine test_sendrecv_2D(npes,pe,root,out_unit)
    integer, intent(in) :: npes,pe,root,out_unit

    if(npes < 3)then
       call mpp_error(FATAL, "Test_sendrecv_2D: minimum of 3 ranks required. Not testing gather; too few ranks.")
    endif
    write(out_unit,*)

    call test_sendrecv_2D_R4(npes, pe, root, out_unit)

    call test_sendrecv_2D_R8(npes, pe, root, out_unit)

    call test_sendrecv_2D_I4(npes, pe, root, out_unit)

    call test_sendrecv_2D_I8(npes, pe, root, out_unit)

  end subroutine test_sendrecv_2D


!> @brief Call the type-specific test_sendrecv_3D routines.
  subroutine test_sendrecv_3D(npes,pe,root,out_unit)
    integer, intent(in) :: npes,pe,root,out_unit

    if(npes < 3)then
       call mpp_error(FATAL, "Test_sendrecv_3D: minimum of 3 ranks required. Not testing gather; too few ranks.")
    endif
    write(out_unit,*)

    call test_sendrecv_3D_R4(npes, pe, root, out_unit)

    call test_sendrecv_3D_R8(npes, pe, root, out_unit)

    call test_sendrecv_3D_I4(npes, pe, root, out_unit)

    call test_sendrecv_3D_I8(npes, pe, root, out_unit)

  end subroutine test_sendrecv_3D


  !> @brief Test together the 2D mpp_send and mpp_recv functions with 32-bit real data arguments.
  subroutine test_sendrecv_2D_R4(npes,pe,root,out_unit)
    integer, intent(in) :: npes,pe,root,out_unit

    integer :: pelist(npes)
    integer :: i,j, p
    real(kind=r4_kind), allocatable, dimension(:,:)  ::  data     !!Data to be sendrecved
    integer :: DS

    DS = 9
    allocate(data(DS, DS))

    !!The full PE list [0, ...,npes-1]
    do i=0,npes-1
       pelist(i+1) = i
    enddo

    !!Initialize all data on all PEs
    data = -1.0
    !! Re-initialize data on the root PE only.
    !! Data is such that we can calculate what it should be with a Formula
    !! using the indecies. E.g.. data(3,4) is 34.000, etc.
    if (pe == root) then
       do i = 1,DS
          do j = 1,DS
             data(i,j) = i*10.0 + j*1.0
          enddo
       enddo
    endif

    !! Send from the source pe all of the data to all other pes
    !! And receive from all other pes
    if ( pe == root ) then
       do p = 1,npes-1
          call mpp_send( data, DS* DS, p )
       end do
    else
       call mpp_recv( data, DS * DS, 0 )
    end if

    call mpp_sync() ! Needed ?

    !! Verify that the data was correctly transmitted.
    if(ANY(pe == pelist(1:npes-1))) then
       do j = 1, DS
          do i = 1, DS
             if (data(i,j) /= ( i*10.0 + j*1.0)) then
                call mpp_error(FATAL, "Test sendrecv 2D R4 failed - basic copy area.")
             endif
          enddo
       enddo
    endif

    call mpp_sync() !
    write(out_unit,*) "Test test_sendrecv_2D_R4  successful ."

  end subroutine test_sendrecv_2D_R4

  !> @brief Test together the 2D mpp_send and mpp_recv functions with 64-bit real data arguments.
  subroutine test_sendrecv_2D_R8(npes,pe,root,out_unit)
      integer, intent(in) :: npes,pe,root,out_unit

    integer :: pelist(npes)
    integer :: i,j, p
    real(kind=r8_kind), allocatable, dimension(:,:)  ::  data     !!Data to be sendrecved
    integer :: DS

    DS = 9
    allocate(data(DS, DS))

    !!The full PE list [0, ...,npes-1]
    do i=0,npes-1
       pelist(i+1) = i
    enddo

    !!Initialize all data on all PEs
    data = -1.0
    !! Re-initialize data on the root PE only.
    !! Data is such that we can calculate what it should be with a Formula
    !! using the indecies. E.g.. data(3,4) is 34.000, etc.
    if (pe == root) then
       do i = 1,DS
          do j = 1,DS
             data(i,j) = i*10.0 + j*1.0
          enddo
       enddo
    endif

    !! Send from the source pe all of the data to all other pes
    !! And receive from all other pes
    if ( pe == root ) then
       do p = 1,npes-1
          call mpp_send( data, DS* DS, p )
       end do
    else
       call mpp_recv( data, DS * DS, 0 )
    end if

    call mpp_sync() ! Needed ?


    !! Verify that the data was correctly transmitted.
    if(ANY(pe == pelist(1:npes-1))) then
       do j = 1, DS
          do i = 1, DS
             if (data(i,j) /= ( i*10.0 + j*1.0)) then
                call mpp_error(FATAL, "Test sendrecv 2D R8 failed - basic copy area.")
             endif
          enddo
       enddo
    endif

    call mpp_sync() !
    write(out_unit,*) "Test test_sendrecv_2D_R8  successful ."

  end subroutine test_sendrecv_2D_R8

  !> @brief Test together the mpp_send and mpp_recv 3D functions with 32-bit real data arguments.
  subroutine test_sendrecv_3D_R4(npes,pe,root,out_unit)
    integer, intent(in) :: npes,pe,root,out_unit

    integer :: pelist(npes)
    integer :: i,j,k, p
    real(kind=r4_kind), allocatable, dimension(:,:,:)  ::  data     !!Data to be sendrecved
    integer :: DS
    integer :: iz, jz  !!The zeroth element to be sendrecved is at pos data(is+iz, js+jz)
    integer :: is, ie, js, je !!The amount of data to be sendrecved is (ie - is)*(je - js)
    integer :: NZ


    NZ = 9 !! Depth of the square tube to be sendrecved.
    DS = 8 !! DS should be less than 10 for the tests below to make sense.
    allocate(data(DS, DS, NZ))

    !!The full PE list is [0, ...,npes-1]
    do i=0,npes-1
       pelist(i+1) = i
    enddo

    !!Initialize all data on all PEs
    data = -1.0
    !! Re-initialize data  on the root PE only.
    !! Data is such that we can calculate what it should be with a Formula
    !! using the indecies. E.g.. data(3,4,5) is 543.000, etc.
    if (pe == root) then
       do i = 1,DS
          do j = 1,DS
             do k = 1,NZ
                data(i,j, k) = k*100.0 + j*10.0 + i*1.0
             enddo
          enddo
       enddo
    endif


    !! Send from the source pe all of the data to all other pes
    !! And receive from all other pes
    if ( pe == root ) then
       do p = 1,npes-1
          call mpp_send( data, DS* DS* NZ, p )
       end do
    else
       call mpp_recv( data, DS * DS * NZ, 0 )
    end if

    call mpp_sync() ! Needed ?

    !! Verify the transmitted data
    if(ANY(pe == pelist(1:npes-1))) then
       !!Note below row (id index of "data() equivalent or formula") changing fastest.
       do k = 1,  NZ
          do j = 1, DS
             do i = 1, DS
                if (data(i,j, k) /= ( k*100.0 + j*10.0 + i*1.0 )) then
                   call mpp_error(FATAL, "Test sendrecv 3D R4 failed - basic copy area.")
                endif
             enddo
          enddo
       enddo
    endif

    call mpp_sync()

    write(out_unit,*) "Test sendrecv 3D R4 successful."

  end subroutine test_sendrecv_3D_R4


  !> @brief Test together the 3D mpp_send and mpp_recv 3D functions with 64-bit real data arguments.
  subroutine test_sendrecv_3D_R8(npes,pe,root,out_unit)
        integer, intent(in) :: npes,pe,root,out_unit

    integer :: pelist(npes)
    integer :: i,j,k, p
    real(kind=r8_kind), allocatable, dimension(:,:,:)  ::  data     !!Data to be sendrecved
    integer :: DS
    integer :: iz, jz  !!The zeroth element to be sendrecved is at pos data(is+iz, js+jz)
    integer :: is, ie, js, je !!The amount of data to be sendrecved is (ie - is)*(je - js)
    integer :: NZ


    NZ = 9 !! Depth of the square tube to be sendrecved.
    DS = 8 !! DS should be less than 10 for the tests below to make sense.
    allocate(data(DS, DS, NZ))

    !!The full PE list is [0, ...,npes-1]
    do i=0,npes-1
       pelist(i+1) = i
    enddo

    !!Initialize all data on all PEs
    data = -1.0
    !! Re-initialize data  on the root PE only.
    !! Data is such that we can calculate what it should be with a Formula
    !! using the indecies. E.g.. data(3,4,5) is 543.000, etc.
    if (pe == root) then
       do i = 1,DS
          do j = 1,DS
             do k = 1,NZ
                data(i,j, k) = k*100.0 + j*10.0 + i*1.0
             enddo
          enddo
       enddo
    endif


    !! Send from the source pe all of the data to all other pes
    !! And receive from all other pes
    if ( pe == root ) then
       do p = 1,npes-1
          call mpp_send( data, DS* DS* NZ, p )
       end do
    else
       call mpp_recv( data, DS * DS * NZ, 0 )
    end if

    call mpp_sync() ! Needed ?

    !! Verify the transmitted data
    if(ANY(pe == pelist(1:npes-1))) then
       !!Note below row (id index of "data() equivalent or formula") changing fastest.
       do k = 1,  NZ
          do j = 1, DS
             do i = 1, DS
                if (data(i,j, k) /= ( k*100.0 + j*10.0 + i*1.0 )) then
                   call mpp_error(FATAL, "Test sendrecv 3D R8 failed - basic copy area.")
                endif
             enddo
          enddo
       enddo
    endif

    call mpp_sync() !

    write(out_unit,*) "Test sendrecv 3D R8 successful."

  end subroutine test_sendrecv_3D_R8

  !> @brief Test together the 2D mpp_send and mpp_recv functions with 32-bit integer data arguments.
  subroutine test_sendrecv_2D_I4(npes,pe,root,out_unit)
    integer, intent(in) :: npes,pe,root,out_unit

    integer :: pelist(npes)
    integer(kind=i4_kind) :: i,j
    integer(kind=i4_kind), allocatable, dimension(:,:)  ::  data     !!Data to be sendrecved
    integer :: DS, p

    DS = 9
    allocate(data(DS, DS))

    !!The full PE list [0, ...,npes-1]
    do i=0,npes-1
       pelist(i+1) = i
    enddo

    !!Initialize all data on all PEs
    data = -1
    !! Re-initialize data on the root PE only.
    !! Data is such that we can calculate what it should be with a Formula
    !! using the indecies. E.g.. data(3,4) is 34.000, etc.
    if (pe == root) then
       do i = 1,DS
          do j = 1,DS
             data(i,j) = i*10 + j
          enddo
       enddo
    endif

    !! Send from the source pe all of the data to all other pes
    !! And receive from all other pes
    if ( pe == root ) then
       do p = 1,npes-1
          call mpp_send( data, DS* DS, p )
       end do
    else
       call mpp_recv( data, DS * DS, 0 )
    end if

    call mpp_sync() ! Needed ?

    !! Verify that the data was correctly transmitted.
    if(ANY(pe == pelist(1:npes-1))) then
       do j = 1, DS
          do i = 1, DS
             if (data(i,j) /= ( i * 10 + j  )) then
                call mpp_error(FATAL, "Test sendrecv 2D I4 failed - basic copy area.")
             endif
          enddo
       enddo
    endif

    call mpp_sync() !
    write(out_unit,*) "Test test_sendrecv_2D_I4  successful ."

  end subroutine test_sendrecv_2D_I4

  !> @brief Test together the 2D mpp_send and mpp_recv functions with 64-bit integer data arguments.
  subroutine test_sendrecv_2D_I8(npes,pe,root,out_unit)
      integer, intent(in) :: npes,pe,root,out_unit

    integer :: pelist(npes)
    integer(kind=i8_kind) :: i,j
    integer(kind=i8_kind), allocatable, dimension(:,:)  ::  data     !!Data to be sendrecved
    integer :: DS, p

    DS = 9
    allocate(data(DS, DS))

    !!The full PE list [0, ...,npes-1]
    do i=0,npes-1
       pelist(i+1) = i
    enddo

    !!Initialize all data on all PEs
    data = -1
    !! Re-initialize data on the root PE only.
    !! Data is such that we can calculate what it should be with a Formula
    !! using the indecies. E.g.. data(3,4) is 34.000, etc.
    if (pe == root) then
       do i = 1,DS
          do j = 1,DS
             data(i,j) = i*10 + j
          enddo
       enddo
    endif

    !! Send from the source pe all of the data to all other pes
    !! And receive from all other pes
    if ( pe == root ) then
       do p = 1,npes-1
          call mpp_send( data, DS* DS, p )
       end do
    else
       call mpp_recv( data, DS * DS, 0 )
    end if

    call mpp_sync() ! Needed ?

    !! Verify that the data was correctly transmitted.
    if(ANY(pe == pelist(1:npes-1))) then
       do j = 1, DS
          do i = 1, DS
             if (data(i,j) /= ( i * 10 + j  )) then
                call mpp_error(FATAL, "Test sendrecv 2D I8 failed - basic copy area.")
             endif
          enddo
       enddo
    endif

    call mpp_sync() !
    write(out_unit,*) "Test test_sendrecv_2D_I8  successful ."

  end subroutine test_sendrecv_2D_I8

  !> @brief Test together the mpp_send and mpp_recv 3D functions with 32-bit integer data arguments.
  subroutine test_sendrecv_3D_I4(npes,pe,root,out_unit)
    integer, intent(in) :: npes,pe,root,out_unit

    integer :: pelist(npes)
    integer(kind=i4_kind) :: i,j,k
    integer(kind=i4_kind), allocatable, dimension(:,:,:)  ::  data     !!Data to be sendrecved
    integer :: DS
    integer :: iz, jz  !!The zeroth element to be sendrecved is at pos data(is+iz, js+jz)
    integer :: is, ie, js, je !!The amount of data to be sendrecved is (ie - is)*(je - js)
    integer :: NZ, p


    NZ = 9 !! Depth of the square tube to be sendrecved.
    DS = 8 !! DS should be less than 10 for the tests below to make sense.
    allocate(data(DS, DS, NZ))

    !!The full PE list is [0, ...,npes-1]
    do i=0,npes-1
       pelist(i+1) = i
    enddo

    !!Initialize all data on all PEs
    data = -1
    !! Re-initialize data  on the root PE only.
    !! Data is such that we can calculate what it should be with a Formula
    !! using the indecies. E.g.. data(3,4,5) is 543.000, etc.
    if (pe == root) then
       do i = 1,DS
          do j = 1,DS
             do k = 1,NZ
                data(i,j, k) = k*100 + j*10 + i
             enddo
          enddo
       enddo
    endif


    !! Send from the source pe all of the data to all other pes
    !! And receive from all other pes
    if ( pe == root ) then
       do p = 1,npes-1
          call mpp_send( data, DS* DS* NZ, p )
       end do
    else
       call mpp_recv( data, DS * DS * NZ, 0 )
    end if

    call mpp_sync() ! Needed ?

    !! Verify the transmitted data
    if(ANY(pe == pelist(1:npes-1))) then
       !!Note below row (id index of "data() equivalent or formula") changing fastest.
       do k = 1,  NZ
          do j = 1, DS
             do i = 1, DS
                if (data(i,j, k) /= ( k * 100 + j*10 + i )) then
                   call mpp_error(FATAL, "Test sendrecv 3D I4 failed - basic copy area.")
                endif
             enddo
          enddo
       enddo
    endif

    call mpp_sync() !

    write(out_unit,*) "Test sendrecv 3D I4 successful."

  end subroutine test_sendrecv_3D_I4

  !> @brief Test together the 3D mpp_send and mpp_recv 3D functions with 64-bit integer data arguments.
  subroutine test_sendrecv_3D_I8(npes,pe,root,out_unit)
        integer, intent(in) :: npes,pe,root,out_unit

    integer :: pelist(npes)
    integer(kind=i8_kind) :: i,j,k
    integer(kind=i8_kind), allocatable, dimension(:,:,:)  ::  data     !!Data to be sendrecved
    integer :: DS
    integer :: iz, jz  !!The zeroth element to be sendrecved is at pos data(is+iz, js+jz)
    integer :: is, ie, js, je !!The amount of data to be sendrecved is (ie - is)*(je - js)
    integer :: NZ, p


    NZ = 9 !! Depth of the square tube to be sendrecved.
    DS = 8 !! DS should be less than 10 for the tests below to make sense.
    allocate(data(DS, DS, NZ))

    !!The full PE list is [0, ...,npes-1]
    do i=0,npes-1
       pelist(i+1) = i
    enddo

    !!Initialize all data on all PEs
    data = -1
    !! Re-initialize data  on the root PE only.
    !! Data is such that we can calculate what it should be with a Formula
    !! using the indecies. E.g.. data(3,4,5) is 543.000, etc.
    if (pe == root) then
       do i = 1,DS
          do j = 1,DS
             do k = 1,NZ
                data(i,j, k) = k*100 + j*10 + i
             enddo
          enddo
       enddo
    endif


    !! Send from the source pe all of the data to all other pes
    !! And receive from all other pes
    if ( pe == root ) then
       do p = 1,npes-1
          call mpp_send( data, DS* DS* NZ, p )
       end do
    else
       call mpp_recv( data, DS * DS * NZ, 0 )
    end if

    call mpp_sync() ! Needed ?

    !! Verify the transmitted data
    if(ANY(pe == pelist(1:npes-1))) then
       !!Note below row (id index of "data() equivalent or formula") changing fastest.
       do k = 1,  NZ
          do j = 1, DS
             do i = 1, DS
                if (data(i,j, k) /= ( k * 100 + j*10 + i )) then
                   call mpp_error(FATAL, "Test sendrecv 3D I8 failed - basic copy area.")
                endif
             enddo
          enddo
       enddo
    endif

    call mpp_sync() !

    write(out_unit,*) "Test sendrecv 3D I8 successful."

  end subroutine test_sendrecv_3D_I8

end program test_mpp_sendrecv
