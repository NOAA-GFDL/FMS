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
!> @brief Test various mpp_gather and mpp_routines.
!> @note  Some of the tested mpp_gather routines are legavy routines originally in file test_mpp.F90.
program test_mpp_gatscat

#ifdef sgi_mipspro
  use shmem_interface
#endif

  use mpp_mod, only : mpp_init, mpp_exit, mpp_pe, mpp_npes, mpp_root_pe, stdout
  use mpp_mod, only : mpp_sync
  use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist, mpp_set_stack_size
  use mpp_mod, only : mpp_gather, mpp_scatter, mpp_error, FATAL
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
  call mpp_set_stack_size(3145746)
  pe = mpp_pe()
  npes = mpp_npes()
  root = mpp_root_pe()
  out_unit = stdout()

  if( pe.EQ.root ) print *, '------------------> Calling test_scatter <------------------'
  call test_scatter_2D(npes,pe,root,out_unit)
  call test_scatter_3D(npes,pe,root,out_unit)
  if( pe.EQ.root ) print *, '------------------> Finished test_scatter <------------------'

  if( pe.EQ.root ) print *, '------------------> Calling test_gather <------------------'
  call test_gather(npes,pe,root,out_unit)
  call test_gatherV(npes,pe,root,out_unit)

  if( pe.EQ.root ) print *, '------------------> Finished test_gather <------------------'

  call MPI_finalize(ierr)


contains

!> @brief    Call some of the type specific (Float vs double) test_scatter_2D routines.
  subroutine test_scatter_2D(npes,pe,root,out_unit)
    integer, intent(in) :: npes,pe,root,out_unit

    if(npes < 3)then
       call mpp_error(FATAL, "Test_scatter_2D: minimum of 3 ranks required. Not testing gather; too few ranks.")
    endif
    write(out_unit,*)

    call test_scatter_2D_R4(npes, pe, root, out_unit)

    call test_scatter_2D_R8(npes, pe, root, out_unit)

  end subroutine test_scatter_2D


!> @brief Call some of the type specific (Float vs double) test_scatter_3D routines.
  subroutine test_scatter_3D(npes,pe,root,out_unit)
    integer, intent(in) :: npes,pe,root,out_unit

    if(npes < 3)then
       call mpp_error(FATAL, "Test_scatter_3D: minimum of 3 ranks required. Not testing gather; too few ranks.")
    endif
    write(out_unit,*)

    call test_scatter_3D_R4(npes, pe, root, out_unit)

    call test_scatter_3D_R8(npes, pe, root, out_unit)

  end subroutine test_scatter_3D


  !> @brief Test the mpp_scatter functions with FLOAT_KIND data arguments.
  subroutine test_scatter_2D_R4(npes,pe,root,out_unit)
    integer, intent(in) :: npes,pe,root,out_unit

    integer :: pelist(npes)
    integer :: i,j,k
    real(kind=r4_kind), allocatable, dimension(:,:)  ::  scatter_data     !!Data to be scattered
    real(kind=r4_kind), allocatable, dimension(:,:)  ::  segment
    integer :: DS, SS  !!Source data size and segment size
    integer :: iz, jz  !!The zeroth element to be scattered is at pos data(is+iz, js+jz)
    integer :: is, ie, js, je !!The amount of data to be scattered is (ie - is)*(je - js)
    integer :: id, jd

    DS = 7 !! DS should be less than 10 for the tests below to make sense.
    SS = 6
    allocate(scatter_data(DS, DS))
    allocate(segment(SS, SS))

    !!The full PE list [0, ...,npes-1]
    do i=0,npes-1
       pelist(i+1) = i
    enddo

    !!Initialize all data on all PEs
    scatter_data = -1
    segment = -2.0
    !! Re-initialize data  on the root PE only.
    !! Data is such that we can calculate what it should be with a Formula
    !! using the indecies. E.g.. data(3,4) is 34.000, etc.
    if (pe == root) then
       do i = 1,DS
          do j = 1,DS
             scatter_data(i,j) = i*10 + j
          enddo
       enddo
       !! And re-initalize segment on the root pe.
       do i = 1,SS
          do j = 1,SS
             segment(i,j) = i * 10 + j
          enddo
       enddo
    endif

    !! Scatter from the source pe a subset of the data array.
    !! The subset is to go into the segment array of the target pes.
    !! The data to scatter is "moved" in a 1D array of size
    !! S=(ie - is) * (je - js) and starts with the data at
    !! position (iz,jz). Recall Fortran is column-major order.
    iz = 2
    jz = 3
    is = 2
    ie = 3
    js = 2
    je = 3
    if(pe .eq. root) then
       call mpp_scatter(is, ie, js, je, pelist(1:npes-1), segment, scatter_data, .true., iz, jz)
    else
       call mpp_scatter(is, ie, js, je, pelist(1:npes -1), segment, scatter_data, .false., iz, jz)
    endif

    call mpp_sync() !


    !! Verify that the segment array has been updated on the target pes (i,e, those
    !! in the pelist, which does not include pe numbered npes)
    if(ANY(pe == pelist(1:npes-1))) then
       i = 1
       j = 1
       !!Note below row (id index of "data() equivalent or formula") changing fastest.
       do jd = js + jz, je + jz
          do id = is + iz, ie + iz
             if (segment(i,j) /= ( id * 10 + jd )) then
                !!write(6,*) i, j, id, jd
                call mpp_error(FATAL, "Test scatter 2D R4 failed in general scatter section.")
             endif
             !! Do to the next data element in segment
             !! If just done the bottom element of a column:
             if(i == SS) then
                i = is
                j = MOD(j + 1, SS) ! next column of segement()
             else
                i = i + 1 ! next row of segemnt()
             endif
          enddo
       enddo
    endif

    call mpp_sync() !
    write(out_unit,*) "Test test_scatter_2D_R4  successful at general scatter section."

    !!Verify that the last pe (numbered npes) did not get the segment array updated!
    if(pe == pelist(npes)) then
       do i = 1,SS
          do j = 1,SS
             if (segment(i,j) /= -2 ) then
                call mpp_error(FATAL, "Test scatter 2D failed. pe=npes segment was changed")
             endif
          end do
       end do
    endif

    call mpp_sync() !
    write(out_unit,*) "Test test_scatter_2D_R4  successful ."

end subroutine test_scatter_2D_R4

  !> @brief Test the mpp_scatter functions with DOUBLE_KIND data arguments.
  subroutine test_scatter_2D_R8(npes,pe,root,out_unit)
    integer, intent(in) :: npes,pe,root,out_unit

    integer :: pelist(npes)
    integer :: i,j,k
    real(kind=r8_kind), allocatable, dimension(:,:)  ::  scatter_data     !!Data to be scattered
    real(kind=r8_kind), allocatable, dimension(:,:)  ::  segment
    integer :: DS, SS  !!Source data size and segment size
    integer :: iz, jz  !!The zeroth element to be scattered is at pos data(is+iz, js+jz)
    integer :: is, ie, js, je !!The amount of data to be scattered is (ie - is)*(je - js)
    integer :: id, jd


    DS = 7 !! DS should be less than 10 for the tests below to make sense.
    SS = 6
    allocate(scatter_data(DS, DS))
    allocate(segment(SS, SS))

    !!The full PE list [0, ...,npes-1]
    do i=0,npes-1
       pelist(i+1) = i
    enddo

    !!Initialize all data on all PEs
    scatter_data = -1
    segment = -2.0
    !! Re-initialize data  on the root PE only.
    !! Data is such that we can calculate what it should be with a Formula
    !! using the indecies. E.g.. data(3,4) is 34.000, etc.
    if (pe == root) then
       do i = 1,DS
          do j = 1,DS
             scatter_data(i,j) = i*10 + j
          enddo
       enddo
       !! And re-initalize segment on the root pe.
       do i = 1,SS
          do j = 1,SS
             segment(i,j) = i * 10 + j
          enddo
       enddo
    endif

    !! Scatter from the source pe a subset of the data array.
    !! The subset is to go into the segment array of the target pes.
    !! The data to scatter is "moved" in a 1D array of size
    !! S=(ie - is) * (je - js) and starts with the data at
    !! position (iz,jz). Recall Fortran is column-major order.
    iz = 2
    jz = 3
    is = 2
    ie = 3
    js = 2
    je = 3
    if(pe .eq. root) then
       call mpp_scatter(is, ie, js, je, pelist(1:npes-1), segment, scatter_data, .true., iz, jz)
    else
       call mpp_scatter(is, ie, js, je, pelist(1:npes -1), segment, scatter_data, .false., iz, jz)
    endif

    call mpp_sync()


    !! Verify that the segment array has been updated on the target pes (i,e, those
    !! in the pelist, which does not include pe numbered npes)
    if(ANY(pe == pelist(1:npes-1))) then
       i = 1
       j = 1
       !!Note below row (id index of "data() equivalent or formula") changing fastest.
       do jd = js + jz, je + jz
          do id = is + iz, ie + iz
             if (segment(i,j) /= ( id * 10 + jd )) then
                !!write(6,*) i, j, id, jd
                call mpp_error(FATAL, "Test scatter 2D R8 failed in general scatter section.")
             endif
             !! Do to the next data element in segment
             !! If just done the bottom element of a column:
             if(i == SS) then
                i = is
                j = MOD(j + 1, SS) ! next column of segement()
             else
                i = i + 1 ! next row of segemnt()
             endif
          enddo
       enddo
    endif

    call mpp_sync()
    write(out_unit,*) "Test test_scatter_2D_R8  successful at general scatter section."

    !!Verify that the last pe (numbered npes) did not get the segment array updated!
    if(pe == pelist(npes)) then
       do i = 1,SS
          do j = 1,SS
             if (segment(i,j) /= -2 ) then
                call mpp_error(FATAL, "Test scatter 2D R8failed. pe=npes segment was changed")
             endif
          end do
       end do
    endif

    call mpp_sync()
    write(out_unit,*) "Test test_scatter_2D_R8  successful ."

end subroutine test_scatter_2D_R8

!> @brief Test the mpp_scatter 3D functions with FLOAT_KIND data arguments.
  subroutine test_scatter_3D_R4(npes,pe,root,out_unit)
    integer, intent(in) :: npes,pe,root,out_unit

    integer :: pelist(npes)
    integer :: i,j,k
    real(kind=r4_kind), allocatable, dimension(:,:,:)  ::  scatter_data     !!Data to be scattered
    real(kind=r4_kind), allocatable, dimension(:,:,:)  ::  segment
    integer :: DS, SS  !!Source data size and segment size
    integer :: iz, jz  !!The zeroth element to be scattered is at pos data(is+iz, js+jz)
    integer :: is, ie, js, je !!The amount of data to be scattered is (ie - is)*(je - js)
    integer :: id, jd, kd
    integer :: NZ
    integer :: dAmount, dCount

    NZ = 11 !! Depth of the square tube to be scattered.
    DS = 6 !! DS should be less than 10 for the tests below to make sense.
    SS = 5 !! Can be different that DS, but see retrictions.
    allocate(scatter_data(DS, DS, NZ))
    allocate(segment(SS, SS, NZ))

    !!The full PE list is [0, ...,npes-1]
    do i=0,npes-1
       pelist(i+1) = i
    enddo

    !!Initialize all data on all PEs
    scatter_data = -1
    segment = -2.0
    !! Re-initialize data  on the root PE only.
    !! Data is such that we can calculate what it should be with a Formula
    !! using the indecies. E.g.. data(3,4,5) is 543.000, etc.
    if (pe == root) then
       do i = 1,DS
          do j = 1,DS
             do k = 1,NZ
             scatter_data(i,j, k) = k*100 + j*10 + i
          enddo
       enddo
       enddo
       !! And re-initalize segment on the root pe.
       do i = 1,SS
          do j = 1,SS
             do k = 1,NZ
             segment(i,j, k) = scatter_data(i,j, k)
          enddo
       enddo
       enddo
    endif

    !! Scatter from the source pe a subset of the data array.
    !! The subset is to go into the segment array of the target pes.
    !! The data to scatter is "moved" in a 1D array of size
    !! S=((ie - is +1) * (je - js + 1) * NZ )and starts with the data at
    !! position (iz,jz, kz). Recall Fortran is column-major order.
    iz = 2
    jz = 2
    is = 2
    ie = 3
    js = 2
    je = 3
    if(pe .eq. root) then
       call mpp_scatter(is, ie, js, je, NZ, pelist(1:npes-1), segment, scatter_data, .true., iz, jz)
    else
       call mpp_scatter(is, ie, js, je, NZ, pelist(1:npes -1), segment, scatter_data, .false., iz, jz)
    endif

    call mpp_sync()

    !! Verify that the segment array has been updated on the target pes (i,e, those
    !! in the pelist, which does not include pe numbered npes)
    !! dAmount is the number of data elements that should be scattered.
    dAmount = (ie - is + 1)*(je -js + 1)*NZ
    dCount = 0;

    if(ANY(pe == pelist(1:npes-1))) then
       kd = 1
       jd = js + jz !(4,5)
       id = is + iz !!increases fastest (4,5)
       !!Note below row (id index of "data() equivalent or formula") changing fastest.
       do k = 1,  NZ
          do j = 1, SS !(je -js + 1)
             do i = 1, SS !(ie - is + 1)
                if(dCount < dAmount) then
                   dCount = dCount + 1
                   !!write(6,*) k, j, i, kd, jd, id
                   if (segment(i,j, k) /= ( kd * 100 + jd*10 + id )) then
                      call mpp_error(FATAL, "Test scatter 3D R4 failed - basic copy area.")
                   endif
                   !! Do to the next data element in segment
                   !!IF the previous one was the corner of a square:
                   if((id == ie + iz  ) .AND.  (jd == (je + jz))) then
                      id = is + iz
                      jd = js + jz
                      kd = kd + 1
                      !!IF the previous one was the botton of a column
                   else if( id == ie + iz  ) then
                      id = is + iz
                      jd = jd + 1
                   else
                      id = id + 1 ! next row of segemnt()
                   endif
                endif
             enddo
          enddo
       enddo
    endif

    call mpp_sync() !
    write(out_unit,*) "Test test_scatter_2D_R4  successful at general scatter section."

    !!Verify that the last pe (numbered npes) did not get the segment array updated!
    if(pe == pelist(npes)) then
       do i = 1,SS
          do j = 1,SS
             do k = 1,  NZ
             if (segment(i,j, k) /= -2 ) then
                call mpp_error(FATAL, "Test scatter 3D R4 failed. pe=npes segment was changed")
             endif
          end do
       end do
       enddo
    endif

    call mpp_sync()
    write(out_unit,*) "Test scatter 3D R4 successful."

  end subroutine test_scatter_3D_R4


  !> @brief Test the mpp_scatter 3D functions with DOUBLE_KIND data arguments.
  subroutine test_scatter_3D_R8(npes,pe,root,out_unit)
    integer, intent(in) :: npes,pe,root,out_unit

    integer :: pelist(npes)
    integer :: i,j,k
    real(kind=r8_kind), allocatable, dimension(:,:,:)  ::  scatter_data     !!Data to be scattered
    real(kind=r8_kind), allocatable, dimension(:,:,:)  ::  segment
    integer :: DS, SS  !!Source data size and segment size
    integer :: iz, jz  !!The zeroth element to be scattered is at pos data(is+iz, js+jz)
    integer :: is, ie, js, je !!The amount of data to be scattered is (ie - is)*(je - js)
    integer :: id, jd, kd
    integer :: NZ
    integer :: dAmount, dCount

    NZ = 11 !! Depth of the square tube to be scattered.
    DS = 6 !! DS should be less than 10 for the tests below to make sense.
    SS = 5 !! Can be different that DS, but see retrictions.
    allocate(scatter_data(DS, DS, NZ))
    allocate(segment(SS, SS, NZ))

    !!The full PE list is [0, ...,npes-1]
    do i=0,npes-1
       pelist(i+1) = i
    enddo

    !!Initialize all data on all PEs
    scatter_data = -1
    segment = -2.0
    !! Re-initialize data  on the root PE only.
    !! Data is such that we can calculate what it should be with a Formula
    !! using the indecies. E.g.. data(3,4,5) is 543.000, etc.
    if (pe == root) then
       do i = 1,DS
          do j = 1,DS
             do k = 1,NZ
             scatter_data(i,j, k) = k*100 + j*10 + i
          enddo
       enddo
       enddo
       !! And re-initalize segment on the root pe.
       do i = 1,SS
          do j = 1,SS
             do k = 1,NZ
             segment(i,j, k) = scatter_data(i,j, k)
          enddo
       enddo
       enddo
    endif

    !! Scatter from the source pe a subset of the data array.
    !! The subset is to go into the segment array of the target pes.
    !! The data to scatter is "moved" in a 1D array of size
    !! S=((ie - is +1) * (je - js + 1) * NZ )and starts with the data at
    !! position (iz,jz, kz). Recall Fortran is column-major order.
    iz = 2
    jz = 2
    is = 2
    ie = 3
    js = 2
    je = 3
    if(pe .eq. root) then
       call mpp_scatter(is, ie, js, je, NZ, pelist(1:npes-1), segment, scatter_data, .true., iz, jz)
    else
       call mpp_scatter(is, ie, js, je, NZ, pelist(1:npes -1), segment, scatter_data, .false., iz, jz)
    endif

    call mpp_sync()

    !! Verify that the segment array has been updated on the target pes (i,e, those
    !! in the pelist, which does not include pe numbered npes)
    !! dAmount is the number of data elements that should be scattered.
    dAmount = (ie - is + 1)*(je -js + 1)*NZ
    dCount = 0;

    if(ANY(pe == pelist(1:npes-1))) then
       kd = 1
       jd = js + jz !(4,5)
       id = is + iz !!increases fastest (4,5)
       !!Note below row (id index of "data() equivalent or formula") changing fastest.
       do k = 1,  NZ
          do j = 1, SS !(je -js + 1)
             do i = 1, SS !(ie - is + 1)
                if(dCount < dAmount) then
                   dCount = dCount + 1
                   !!write(6,*) k, j, i, kd, jd, id
                   if (segment(i,j, k) /= ( kd * 100 + jd*10 + id )) then
                      call mpp_error(FATAL, "Test scatter 3D R8 failed - basic copy area.")
                   endif
                   !! Do to the next data element in segment
                   !!IF the previous one was the corner of a square:
                   if((id == ie + iz  ) .AND.  (jd == (je + jz))) then
                      id = is + iz
                      jd = js + jz
                      kd = kd + 1
                      !!IF the previous one was the botton of a column
                   else if( id == ie + iz  ) then
                      id = is + iz
                      jd = jd + 1
                   else
                      id = id + 1 ! next row of segemnt()
                   endif
                endif
             enddo
          enddo
       enddo
    endif

    call mpp_sync() !
    write(out_unit,*) "Test test_scatter_2D_R8  successful at general scatter section."

    !!Verify that the last pe (numbered npes) did not get the segment array updated!
    if(pe == pelist(npes)) then
       do i = 1,SS
          do j = 1,SS
             do k = 1,  NZ
             if (segment(i,j, k) /= -2 ) then
                call mpp_error(FATAL, "Test scatter 3D R8 failed. pe=npes segment was changed")
             endif
          end do
       end do
       enddo
    endif

    call mpp_sync()
    write(out_unit,*) "Test scatter 3D R8 successful."

  end subroutine test_scatter_3D_R8



  !> @brief  Call some of the type specific (Float vs double) test_gather routines.
    subroutine test_gather(npes,pe,root,out_unit)
     integer, intent(in) :: npes,pe,root,out_unit

     if(npes < 3)then
       call mpp_error(FATAL, "Test_gather: minimum of 3 ranks required. Not testing gather; too few ranks.")
     endif
     write(out_unit,*)

     call test_gather_R4(npes, pe, root, out_unit)

     call test_gather_R8(npes, pe, root, out_unit)

   end subroutine test_gather

  !> @brief  Test the scalar mpp_gather routine with FLOAT_KIND data.
     subroutine test_gather_R4(npes,pe,root,out_unit)
     integer, intent(in) :: npes,pe,root,out_unit

     integer :: pelist(npes)
     integer :: i
     real(kind=r4_kind) :: rdata(npes)
     real(kind=r4_kind) :: val

     if(npes < 3)then
       call mpp_error(FATAL, "Test_gather: minimum of 3 ranks required. Not testing gather; too few ranks.")
     endif
     write(out_unit,*)

     val = pe
     rdata = -1.0
     do i=1,npes
       pelist(i) = i-1
     enddo

     call mpp_gather((/val/),rdata)
     if(pe == root)then
       do i=1,npes
        if(INT(rdata(i)) /= pelist(i))then
           write(6,*) "Gathered data ",INT(rdata(i)), " NE reference ",pelist(i), "at i=",i
           call mpp_error(FATAL, "Test gather R4 uniform vector with global pelist failed")
        endif
       enddo
     endif

     call mpp_sync()
     write(out_unit,*) "Test gather uniform vector with global pelist successful"

     rdata = -1.0
     if(ANY(pe == pelist(2:npes)))call mpp_gather((/val/),rdata(2:npes),pelist(2:npes))
     if(pe == pelist(2))then
       do i=2,npes
        if(INT(rdata(i)) /= pelist(i))then
           write(6,*) "Gathered data ",INT(rdata(i)), " NE reference ",pelist(i), "at i=",i
           call mpp_error(FATAL, "Test gather R4 uniform vector with reduced pelist failed")
        endif
       enddo
     endif
     call mpp_sync()
     write(out_unit,*) "Test gather uniform vector with reduced pelist successful"

  end subroutine test_gather_R4


  !> @brief  Test the scalar mpp_gather routine with DOUBLE_KIND data.
    subroutine test_gather_R8(npes,pe,root,out_unit)
     integer, intent(in) :: npes,pe,root,out_unit

     integer :: pelist(npes)
     integer :: i
     real(kind=r8_kind) :: rdata(npes)
     real(kind=r8_kind) :: val

     if(npes < 3)then
       call mpp_error(FATAL, "Test_gather: minimum of 3 ranks required. Not testing gather; too few ranks.")
     endif
     write(out_unit,*)

     val = pe
     rdata = -1.0
     do i=1,npes
       pelist(i) = i-1
     enddo

     call mpp_gather((/val/),rdata)
     if(pe == root)then
       do i=1,npes
        if(INT(rdata(i)) /= pelist(i))then
           write(6,*) "Gathered data ",INT(rdata(i)), " NE reference ",pelist(i), "at i=",i
           call mpp_error(FATAL, "Test gather R8 uniform vector with global pelist failed")
        endif
       enddo
     endif

     call mpp_sync()
     write(out_unit,*) "Test gather uniform vector with global pelist successful"

     rdata = -1.0
     if(ANY(pe == pelist(2:npes)))call mpp_gather((/val/),rdata(2:npes),pelist(2:npes))
     if(pe == pelist(2))then
       do i=2,npes
        if(INT(rdata(i)) /= pelist(i))then
           write(6,*) "Gathered data ",INT(rdata(i)), " NE reference ",pelist(i), "at i=",i
           call mpp_error(FATAL, "Test gather R8 uniform vector with reduced pelist failed")
        endif
       enddo
     endif
     call mpp_sync()
     write(out_unit,*) "Test gather uniform vector with reduced pelist successful"

   end subroutine test_gather_R8

   !> @brief  Test the 1Dvector  mpp_gather routine.
   !> @todo   Change or refactor this routine to explicitly use FLOAT_KIND and DOUBLE_KIND.
  subroutine test_gatherV(npes,pe,root,out_unit)
  implicit none
     integer, intent(in) :: npes,pe,root,out_unit

     integer :: pelist(npes),rsize(npes)
     integer :: i,j,k,dsize,ssize
     real,allocatable :: sdata(:), rdata(:), ref(:)

     if(npes < 3)then
       call mpp_error(FATAL, "Test_gatherV: minimum of 3 ranks required. Not testing gather; too few ranks.")
     elseif(npes > 9999)then
       call mpp_error(FATAL, "Test_gatherV: maximum of 9999 ranks supported. Not testing gatherV; too many ranks.")
     endif
     write(out_unit,*)

     ssize = pe+1
     allocate(sdata(ssize))
     do i=1,ssize
       sdata(i) = pe + 0.0001*i
     enddo
     do i=1,npes
       pelist(i) = i-1
       rsize(i) = i
     enddo

     dsize = sum(rsize)
     allocate(rdata(dsize),ref(dsize))
     rdata = -1.0
     k=1
     do j=1,npes
       do i=1,rsize(j)
          ref(k) = pelist(j) + 0.0001*i
          k = k+1
     enddo;enddo

     call mpp_gather(sdata,ssize,rdata,rsize)

     if(pe == root)then
       k = 1
       do j=1,npes
         do i=1,rsize(j)
           if(rdata(k) /= ref(k))then
              write(6,*) "Gathered data ",rdata(k), " NE reference ",ref(k), "at k=",k
              call mpp_error(FATAL, "Test gatherV global pelist failed")
           endif
           k = k+1
       enddo;enddo
     endif

     call mpp_sync()
     write(out_unit,*) "Test gatherV with global pelist successful"

     rdata = -1.0
     ref(1) = -1.0

     if(ANY(pe == pelist(2:npes)))call mpp_gather(sdata,ssize,rdata(2:),rsize(2:),pelist(2:npes))

     if(pe == pelist(2))then
       k = 1
       do j=1,npes
         do i=1,rsize(j)
           if(rdata(k) /= ref(k) )then
              write(6,*) "Gathered data ",rdata(k), " NE reference ",ref(k), "at k=",k
              call mpp_error(FATAL, "Test gatherV with reduced pelist failed")
           endif
           k = k+1
       enddo;enddo
     endif
     call mpp_sync()

     write(out_unit,*) "Test gatherV with reduced pelist successful"
     deallocate(sdata,rdata,ref)
  end subroutine test_gatherV

end program test_mpp_gatscat
