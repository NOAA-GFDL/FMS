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
!> @author Ryan Mulhall
!> @email gfdl.climate.model.info@noaa.gov
!> @brief Unit tests for mpp global max, min, and sum functions
!> @description Generates a random data set for both SIZEs of reals and ints
!> then checks routines with local results received from each pe
program test_global_arrays

  use platform_mod
  use mpp_mod,         only: mpp_init, mpp_exit, mpp_pe, mpp_npes, mpp_root_pe
  use mpp_mod,         only: mpp_set_stack_size, mpp_sync, mpp_sync_self
  use mpp_mod,         only: mpp_error, FATAL, NOTE, mpp_send, mpp_recv, WARNING
  use mpp_mod,         only: mpp_init_test_init_true_only, mpp_set_root_pe
  use mpp_io_mod,      only: mpp_io_init
  use mpp_domains_mod, only: mpp_domains_init, mpp_define_domains, domain2d
  use mpp_domains_mod, only: mpp_define_layout, mpp_domains_set_stack_size
  use mpp_domains_mod, only: mpp_get_global_domain, mpp_global_max
  use mpp_domains_mod, only: mpp_global_min, mpp_get_data_domain,mpp_get_compute_domain
  use mpp_domains_mod, only: mpp_domains_exit, mpp_update_domains
  use mpp_domains_mod, only: mpp_get_domain_shift, mpp_global_sum

  implicit none

  integer, parameter            :: length=64
  integer                       :: id, pe, npes, root, i, j, icount, jcount
  integer(i4_kind)              :: maxI4, minI4, ierr, sumI4, sumI4_5d
  integer(i8_kind)              :: maxI8, minI8, sumI8, sumI8_5d
  integer(i4_kind), allocatable :: dataI4(:,:), dataI4_5d(:,:,:,:,:), dataI4_shuf(:,:)
  integer(i8_kind), allocatable :: dataI8(:,:), dataI8_5d(:,:,:,:,:), dataI8_shuf(:,:)
  real(r4_kind), allocatable    :: dataR4(:,:), dataR4_5d(:,:,:,:,:), dataR4_shuf(:,:)
  real(r8_kind), allocatable    :: dataR8(:,:), dataR8_5d(:,:,:,:,:), dataR8_shuf(:,:)
  real, allocatable             :: rands(:)
  type(domain2D)                :: domain
  real(r8_kind)                 :: rcoef, maxR8, minR8, sumR8, sumR8_5d
  real(r4_kind)                 :: maxR4, minR4, sumR4, sumR4_5d
  integer                       :: isc, iec, jsc, jec
  integer                       :: isd, ied, jsd, jed
  character(len=32)             :: strTmp1, strTmp2
  integer(i4_kind), parameter   :: randmaxI4 = 2048
  integer(i8_kind), parameter   :: randmaxI8 = 4096

  call mpp_init(mpp_init_test_init_true_only)
  call mpp_io_init()
  call mpp_domains_init()
  call mpp_set_stack_size(3145746)
  call mpp_domains_set_stack_size(3145746)
  pe = mpp_pe()
  npes = mpp_npes()
  call mpp_set_root_pe(0)
  root = mpp_root_pe()

  !> define domains and allocate
  call mpp_define_domains( (/1,length,1,length/), (/4,2/), domain, xhalo=0)
  call mpp_get_compute_domain(domain, jsc, jec, isc, iec)
  call mpp_get_data_domain(domain, jsd, jed, isd, ied)
  allocate(dataI4(jsd:jed, isd:ied),dataI8(jsd:jed, isd:ied), rands(length*length))
  allocate(dataR4(jsd:jed, isd:ied), dataR8(jsd:jed, isd:ied))
  allocate(dataR4_shuf(jsd:jed, isd:ied), dataR8_shuf(jsd:jed, isd:ied))
  allocate(dataI4_shuf(jsd:jed, isd:ied), dataI8_shuf(jsd:jed, isd:ied))

  dataI4 = 0; dataI8 = 0; dataR4 = 0.0; dataR8 = 0.0
  dataR8_shuf=0.0; dataR4_shuf=0.0;dataI8_shuf=0; dataI4_shuf=0

  !> make random arrays
  call random_seed()
  call random_number(rands)
  do i=isc, iec-1
    do j=jsc, jec-1
      rcoef = rands(j + i*length) * 2 -1
      dataI4(j, i) = int(rcoef * randmaxI4, kind=i4_kind)
      dataI8(j, i) = int(rcoef * randmaxI8, kind=i8_kind)
      dataR4(j, i) = real(rcoef, kind=r4_kind)
      dataR8(j, i) = real(rcoef, kind=r8_kind)
    end do
  end do

  !> test global max and mins from each kind
  call mpp_error(NOTE, "----------Testing 32-bit int mpp_global_max and mpp_global_min----------")
  call mpp_update_domains(dataI4, domain)
  maxI4 = mpp_global_max(domain, dataI4)
  minI4 = mpp_global_min(domain, dataI4)
  write(strTmp1, *) maxI4
  write(strTmp2, *) minI4
  if(.NOT. checkResultInt4((/minI4, maxI4 /))) then
    call mpp_error(FATAL, "test_global_arrays: invalid 32-bit integer results"// &
                               NEW_LINE('a')//"Max: "//strTmp1//" Min: "//strTmp2 )
  endif
  call mpp_sync()
  call mpp_error(NOTE, "----------Testing 64-bit int mpp_global_max and mpp_global_min----------")
  call mpp_update_domains(dataI8, domain)
  maxI8 = mpp_global_max(domain, dataI8)
  minI8 = mpp_global_min(domain, dataI8)
  write(strTmp1, *) maxI8
  write(strTmp2, *) minI8
  if(.NOT. checkResultInt8((/minI8, maxI8 /))) then
    call mpp_error(FATAL, "test_global_arrays: invalid 64-bit integer results"// &
                               NEW_LINE('a')//"Max: "//strTmp1//" Min: "//strTmp2 )
  endif
  call mpp_sync()
  call mpp_error(NOTE, "----------Testing 32-bit real mpp_global_max and mpp_global_min----------")
  call mpp_update_domains(dataR4, domain)
  maxR4 = mpp_global_max(domain, dataR4)
  minR4 = mpp_global_min(domain, dataR4)
  write(strTmp1, *) maxR4
  write(strTmp2, *) minR4
  if(.NOT. checkResultReal4((/minR4, maxR4 /))) then
    call mpp_error(FATAL, "test_global_arrays: invalid 32-bit real results"// &
                               NEW_LINE('a')//"Max: "//strTmp1//" Min: "//strTmp2 )
  endif
  call mpp_sync()
  call mpp_error(NOTE, "----------Testing 64-bit real mpp_global_max and mpp_global_min----------")
  call mpp_update_domains(dataR8, domain)
  maxR8 = mpp_global_max(domain, dataR8)
  minR8 = mpp_global_min(domain, dataR8)
  write(strTmp1, *) maxR8
  write(strTmp2, *) minR8
  if(.NOT. checkResultReal8((/minR8, maxR8 /))) then
    call mpp_error(FATAL, "test_global_arrays: invalid 64-bit real results"// &
                               NEW_LINE('a')//"Max: "//strTmp1//" Min: "//strTmp2 )
  endif

  !> test global sums for each kind
  call mpp_error(NOTE, "----------Testing 32-bit real mpp_global_sum----------")
  call mpp_update_domains(dataR4, domain)
  sumR4 = mpp_global_sum(domain, dataR4)
  write(strTmp1,*) sumR4
  if(.NOT. checkSumReal4(sumR4)) then
    call mpp_error(FATAL, "test_global_arrays: invalid 32-bit real sum"// &
                               NEW_LINE('a')//"Sum: "// strTmp1 )
  endif
  call mpp_error(NOTE, "----------Testing 64-bit real mpp_global_sum----------")
  call mpp_update_domains(dataR8, domain)
  sumR8 = mpp_global_sum(domain, dataR8)
  write(strTmp1,*) sumR8
  if(.NOT. checkSumReal8(sumR8)) then
    call mpp_error(FATAL, "test_global_arrays: invalid 64-bit real sum"// &
                               NEW_LINE('a')//"Sum: "// strTmp1 )
  endif
  call mpp_error(NOTE, "----------Testing 32-bit integer mpp_global_sum----------")
  call mpp_update_domains(dataI4, domain)
  sumI4 = mpp_global_sum(domain, dataI4)
  write(strTmp1,*) sumI4
  if(.NOT. checkSumInt4(sumI4)) then
    call mpp_error(FATAL, "test_global_arrays: invalid 32-bit integer sum"// &
                               NEW_LINE('a')//"Sum: "// strTmp1 )
  endif
  call mpp_error(NOTE, "----------Testing 64-bit integer mpp_global_sum----------")
  call mpp_update_domains(dataI8, domain)
  sumI8 = mpp_global_sum(domain, dataI8)
  write(strTmp1,*) sumI8
  if(.NOT. checkSumInt8(sumI8)) then
    call mpp_error(FATAL, "test_global_arrays: invalid 64-bit integer sum"// &
                               NEW_LINE('a')//"Sum: "// strTmp1 )
  endif

  !> shuffle real data ordering and copy into array with 5 ranks
  dataR4_shuf = dataR4
  dataR8_shuf = dataR8
  call shuffleDataR4(dataR4_shuf)
  call shuffleDataR8(dataR8_shuf)
  allocate(dataR4_5d(jsd:jed, isd:ied, 1, 1, 1), dataR8_5d(jsd:jed,isd:ied, 1, 1, 1))

  dataR4_5d = 0.0
  dataR8_5d = 0.0

  do i=isc,iec
    do j=jsc,jec
       dataR4_5d(j, i, 1, 1, 1) = dataR4_shuf(j, i)
       dataR8_5d(j, i, 1, 1, 1) = dataR8_shuf(j, i)
    end do
  end do
  call mpp_sync()

  call mpp_error(NOTE, "----------Testing 32-bit real mpp_global_sum with 5 ranks and reordering----------")
  call mpp_update_domains(dataR4_5d, domain)
  sumR4_5d = mpp_global_sum(domain, dataR4_5d)
  ! check that shuffled array results are approximately the same as the original array
  if(abs(sumR4-sumR4_5d) .gt. 1E-4 ) then
    strTmp1 = ""; strTmp2=""
    write(strTmp1,*) sumR4_5d
    write(strTmp2,*) sumR4
    call mpp_error(FATAL,"test_global_arrays: invalid 32-bit real answer after reordering"// &
                   NEW_LINE('a')//"Sum: "// strTmp1// " ne "//strTmp2)
  endif

  call mpp_error(NOTE, "----------Testing 64-bit real mpp_global_sum with 5 ranks and reordering----------")
  call mpp_update_domains(dataR8_5d, domain)
  sumR8_5d = mpp_global_sum(domain, dataR8_5d)
  ! check that shuffled array results are approximately the same as the original array
  if(abs(sumR8-sumR8_5d) .gt. 1E-7) then
    strTmp1 = ""; strTmp2=""
    write(strTmp1,*) sumR8_5d
    write(strTmp2,*) sumR8
    call mpp_error(FATAL,"test_global_arrays: invalid 64-bit real answer after reordering"// &
                   NEW_LINE('a')//"Sum: "// strTmp1// " ne "//strTmp2)
  endif

  !> shuffle integer data ordering and copy into array with 5 ranks
  dataI4_shuf = dataI4
  dataI8_shuf = dataI8
  call shuffleDataI4(dataI4_shuf)
  call shuffleDataI8(dataI8_shuf)
  allocate(dataI4_5d(jsd:jed, isd:ied, 1, 1, 1), dataI8_5d(jsd:jed,isd:ied, 1, 1, 1))

  dataI4_5d = 0
  dataI8_5d = 0
  do i=isc,iec
    do j=jsc,jec
      dataI4_5d(j, i, 1, 1, 1) = dataI4_shuf(j, i)
      dataI8_5d(j, i, 1, 1, 1) = dataI8_shuf(j, i)
    end do
  end do
  call mpp_sync()

  call mpp_error(NOTE, "----------Testing 32-bit integer mpp_global_sum with 5 ranks and reordering----------")
  call mpp_update_domains(dataI4_5d, domain)
  sumI4_5d = mpp_global_sum(domain, dataI4_5d)

  ! check that shuffled array results are approximately the same as the original array
  if(sumI4 .ne. sumI4_5d) then
    strTmp1 = ""; strTmp2=""
    write(strTmp1,*) sumI4_5d
    write(strTmp2,*) sumI4
    call mpp_error(FATAL,"test_global_arrays: invalid 32-bit integer answer after reordering"// &
                   NEW_LINE('a')//"Sum: "// strTmp1// " ne "//strTmp2)
  endif

  call mpp_error(NOTE, "----------Testing 64-bit integer mpp_global_sum with 5 ranks and reordering----------")
  call mpp_update_domains(dataI8_5d, domain)
  sumI8_5d = mpp_global_sum(domain, dataI8_5d)

  ! check that shuffled array results are approximately the same as the original array
  !> @note This test fails with gcc 9.3.0
  if(sumI8 .ne. sumI8_5d) then
    strTmp1 = ""; strTmp2=""
    write(strTmp1,*) sumI8_5d
    write(strTmp2,*) sumI8
    call mpp_error(FATAL,"test_global_arrays: invalid 64-bit integer answer after reordering"// &
                   NEW_LINE('a')//"Sum: "// strTmp1// " ne "//strTmp2)
  endif

  deallocate(dataI4, dataI8, dataR4, dataR8, rands, dataI4_5d, dataI8_5d, dataR4_5d, dataR8_5d)
  deallocate(dataR4_shuf, dataR8_shuf,dataI4_shuf, dataI8_shuf)
  call mpp_domains_exit()
  call MPI_FINALIZE(ierr)

  contains

!> true if all pes return the same result and have a lower/higher local max/min
function checkResultInt4(res)
  logical                               :: checkResultInt4
  integer(i4_kind),intent(in)           :: res(2)
  integer(i4_kind),allocatable          :: tres(:)

  allocate(tres(2))
  checkResultInt4 = res(2).GE.maxval(dataI4) .and. res(1).LE.minval(dataI4)
  if(.NOT.checkResultInt4) then
    return
  end if
  !> check that all pes have same results
  if( pe.EQ.root) then
    tres = res
    do i=1, npes-1
      call mpp_send(tres,2, i)
    end do
    checkResultInt4 = .true.
  else
    call mpp_recv(tres,2, root)
    checkResultInt4 = checkResultInt4 .and. res(1) .EQ. tres(1) .and. res(2) .eq. tres(2)
  end if
  call mpp_sync()
  deallocate(tres)
end function checkResultInt4

function checkResultInt8(res)
  logical                               :: checkResultInt8
  integer(i8_kind),intent(in)           :: res(2)
  integer(i8_kind),allocatable          :: tres(:)

  allocate(tres(2))
  checkResultInt8 = res(2).GE.maxval(dataI8) .and. res(1).LE.minval(dataI8)
  if(.NOT.checkResultInt8) then
    return
  end if
  !> check that all pes have same results
  if( pe.EQ.root) then
    tres = res
    do i=1, npes-1
      call mpp_send(tres,2, i)
    end do
    checkResultInt8 = .true.
  else
    call mpp_recv(tres,2, root)
    checkResultInt8 = checkResultInt8 .and. res(1) .EQ. tres(1) .and. res(2) .eq. tres(2)
  end if
  call mpp_sync()
  deallocate(tres)
end function checkResultInt8

function checkResultReal4(res)
  logical                            :: checkResultReal4
  real(r4_kind),intent(in)           :: res(2)
  real(r4_kind),allocatable          :: tres(:)

  allocate(tres(2))
  checkResultReal4 = res(2).GE.maxval(dataR4) .and. res(1).LE.minval(dataR4)
  if(.NOT. checkResultReal4) then
    return
  end if
  !> check that all pes have same results
  if( pe.EQ.root) then
    tres = res
    do i=1, npes-1
      call mpp_send(tres,2, i)
    end do
    checkResultReal4 = .true.
  else
    call mpp_recv(tres,2, root)
    checkResultReal4 = checkResultReal4 .and. (abs((res(1)-tres(1))/res(1)) .lt. 1e-4) .and. &
                       (abs((res(2)-tres(2))/res(2)) .lt. 1e-4)
  end if
  call mpp_sync()
  deallocate(tres)
end function checkResultReal4

function checkResultReal8(res)
  logical                            :: checkResultReal8
  real(r8_kind),intent(in)           :: res(:)
  real(r8_kind),allocatable          :: tres(:)

  allocate(tres(2))
  checkResultReal8 = res(2).GE.maxval(dataR8) .and. res(1).LE.minval(dataR8)
  if(.NOT.checkResultReal8) then
    return
  end if
  !> check that all pes have same results
  if( pe.EQ.root) then
    tres = res
    do i=1, npes-1
      call mpp_send(tres,2, i)
    end do
    checkResultReal8 = .true.
  else
    call mpp_recv(tres,2, root)
    checkResultReal8 = checkResultReal8 .and. (abs((res(1)-tres(1))/res(1)) .lt. 1e-7) .and. &
                       (abs((res(2)-tres(2))/res(2)) .lt. 1e-7)
  end if
  call mpp_sync()
  deallocate(tres)
end function checkResultReal8

!>@brief Sum local sums from pes and compares with gsum
!>@return True if gsum is the global sum, false otherwise
function checkSumReal4(gsum)
  logical                   :: checkSumReal4
  real(r4_kind),intent(in)  :: gsum
  real(r4_kind),allocatable :: recv(:) !> pe's local sum at 1, global sum at 2
  real(r4_kind)             :: nsum
  integer                   :: i

  allocate(recv(2))
  ! root receives and sums local sums from each pe
  if(pe .eq. root) then
    nsum = SUM(dataR4)
    do i=1, npes - 1
      call mpp_recv(recv, 2, i)
      nsum = nsum + recv(1)
      ! also check for matching global sum
      if( abs((recv(2)-gsum)/gsum) .gt. 1e-4) then
        checkSumReal4 = .false.
        deallocate(recv)
        return
      endif
    end do
    checkSumReal4 = (abs((nsum-gsum)/gsum) .lt. 1e-4)
  else
    recv(1) = SUM(dataR4)
    recv(2) = gsum
    call mpp_send(recv, 2, root)
    checkSumReal4 = .true.
  endif
  call mpp_sync()
  deallocate(recv)
end function checkSumReal4

!>@brief Sum local sums from pes and compares with gsum
!>@return True if gsum is the global sum, false otherwise
function checkSumReal8(gsum)
  logical                   :: checkSumReal8
  real(r8_kind),intent(in)  :: gsum
  real(r8_kind),allocatable :: recv(:) !> pe's local sum at 1, global sum at 2
  real(r8_kind)             :: nsum
  integer                   :: i

  allocate(recv(2))
  ! root receives and sums local sums from each pe
  if(pe .eq. root) then
    nsum = SUM(dataR8)
    do i=1, npes - 1
      call mpp_recv(recv, 2, i)
      nsum = nsum + recv(1)
      ! also check for matching global sum
      if( abs((recv(2)-gsum)/gsum) .gt. 1e-7 ) then
        checkSumReal8 = .false.
        deallocate(recv)
        return
      endif
    end do
    checkSumReal8 = (abs((nsum-gsum)/gsum) .lt. 1e-7)
  else
    recv(1) = SUM(dataR8)
    recv(2) = gsum
    call mpp_send(recv, 2, root)
    checkSumReal8 = .true.
  endif
  call mpp_sync()
  deallocate(recv)
end function checkSumReal8

!>@brief Sum local sums from pes and compares with gsum
!>@return True if gsum is the global sum, false otherwise
function checkSumInt4(gsum)
  logical                      :: checkSumInt4
  integer(i4_kind),intent(in)  :: gsum
  integer(i4_kind),allocatable :: recv(:) !> pe's local sum at 1, global sum at 2
  integer(i4_kind)             :: nsum
  integer                      :: i

  allocate(recv(2))
  ! root receives and sums local sums from each pe
  if(pe .eq. root) then
    nsum = SUM(dataI4)
    do i=1, npes - 1
      call mpp_recv(recv, 2, i)
      nsum = nsum + recv(1)
      ! also check for matching global sum
      if( recv(2) .ne. gsum ) then
        checkSumInt4 = .false.
        deallocate(recv)
        return
      endif
    end do
    checkSumInt4 = nsum .eq. gsum
  else
    recv(1) = SUM(dataI4)
    recv(2) = gsum
    call mpp_send(recv, 2, root)
    checkSumInt4 = .true.
  endif
  call mpp_sync()
  deallocate(recv)
end function checkSumInt4

!>@brief Sum local sums from pes and compares with gsum
!>@return True if gsum is the global sum, false otherwise
function checkSumInt8(gsum)
  logical                      :: checkSumInt8
  integer(i8_kind),intent(in)  :: gsum
  integer(i8_kind),allocatable :: recv(:) !> pe's local sum at 1, global sum at 2
  integer(i8_kind)             :: nsum
  integer                      :: i

  allocate(recv(2))
  ! root receives and sums local sums from each pe
  if(pe .eq. root) then
    nsum = SUM(dataI8)
    do i=1, npes - 1
      call mpp_recv(recv, 2, i)
      nsum = nsum + recv(1)
      ! also check for matching global sum
      if( recv(2) .ne. gsum ) then
        checkSumInt8 = .false.
        deallocate(recv)
        return
      endif
    end do
    checkSumInt8 = nsum .eq. gsum
  else
    recv(1) = SUM(dataI8)
    recv(2) = gsum
    call mpp_send(recv, 2, root)
    checkSumInt8 = .true.
  endif
  call mpp_sync()
  deallocate(recv)
end function checkSumInt8

!> aggregates data on root and randomizes ordering, then sends partitions back to pes
subroutine shuffleDataI4(dataI4)
  integer(i4_kind), intent(INOUT) :: dataI4(:,:)
  integer(i4_kind), allocatable :: trans(:,:), shuffled(:),tmp
  integer :: rind

  allocate(trans(SIZE(dataI4,1), SIZE(dataI4,2)))
  allocate(shuffled(1:length*length))

  if( pe.eq.root) then
    !> get array partitions and aggregate into 1d
    shuffled(1:SIZE(dataI4)) = RESHAPE(dataI4, (/SIZE(dataI4)/))
    do i=1, npes-1
      call mpp_recv(trans, SIZE(dataI4) , i)
      shuffled( SIZE(trans)*i+1 : SIZE(trans)*(i+1)) = RESHAPE(trans, (/SIZE(trans)/))
    end do

    !> shuffle order
    do i=1, length*length
      rind = (rands(i) * length * length)
      if( rind .eq. 0) then
        rind = 1
      endif
      tmp = shuffled(i)
      shuffled(i) = shuffled(rind)
      shuffled(rind) = tmp
    end do
    trans = 0

    !> send back to pes
    do i=0, npes-1
      trans = RESHAPE(shuffled(SIZE(trans)*i + 1:SIZE(trans)*(i+1)), &
                         (/SIZE(trans,1), SIZE(trans,2) /) )
      if(i.ne.root) then
        call mpp_send(trans, SIZE(trans), i)
      else
        dataI4 = trans
      endif
    end do
  else
    call mpp_send(dataI4, SIZE(dataI4), root)
    call mpp_recv(trans, SIZE(dataI4), root)
    dataI4 = trans
  endif
  deallocate(trans, shuffled)
end subroutine shuffleDataI4

!> aggregates data on root and randomizes ordering, then sends partitions back to pes
subroutine shuffleDataI8(dataI8)
  integer(i8_kind), intent(INOUT) :: dataI8(:,:)
  integer(i8_kind), allocatable :: trans(:,:), shuffled(:), tmp
  integer :: rind

  allocate(trans(SIZE(dataI8,1), SIZE(dataI8,2)))
  allocate(shuffled(1:length*length))

  if( pe.eq.root) then
    !> get array partitions and aggregate into 1d
    shuffled(1:SIZE(dataI8)) = RESHAPE(dataI8, (/SIZE(dataI8)/))
    do i=1, npes-1
      call mpp_recv(trans, SIZE(dataI8) , i)
      shuffled( SIZE(trans)*i+1 : SIZE(trans)*(i+1)) = RESHAPE(trans, (/SIZE(trans)/))
    end do

    !> shuffle order
    do i=1, length*length
      rind = (rands(i) * length * length)
      if( rind .eq. 0) then
        rind = 1
      endif
      tmp = shuffled(i)
      shuffled(i) = shuffled(rind)
      shuffled(rind) = tmp
    end do
    trans = 0

    !> send back to pes
    do i=0, npes-1
      trans = RESHAPE(shuffled(SIZE(trans)*i + 1:SIZE(trans)*(i+1)), &
                         (/SIZE(trans,1), SIZE(trans,2) /) )
      if(i.ne.root) then
        call mpp_send(trans, SIZE(trans), i)
      else
        dataI8 = trans
      endif
    end do
  else
    call mpp_send(dataI8, SIZE(dataI8), root)
    call mpp_recv(trans, SIZE(dataI8), root)
    dataI8 = trans
  endif
  deallocate(trans, shuffled)
end subroutine shuffleDataI8

!> aggregates 32-bit real data on root and randomizes ordering, then sends partitions back to pes
subroutine shuffleDataR4(dataR4)
  real(r4_kind), intent(INOUT) :: dataR4(:,:)
  real(r4_kind), allocatable :: trans(:,:), shuffled(:), tmp
  integer :: rind

  allocate(trans(SIZE(dataR4,1), SIZE(dataR4,2)))
  allocate(shuffled(1:length*length))

  if( pe.eq.root) then
    !> get array partitions and aggregate into 1d
    shuffled(1:SIZE(dataR4)) = RESHAPE(dataR4, (/SIZE(dataR4)/))
    do i=1, npes-1
      call mpp_recv(trans, SIZE(dataR4) , i)
      shuffled( SIZE(trans)*i+1 : SIZE(trans)*(i+1)) = RESHAPE(trans, (/SIZE(trans)/))
    end do

    !> shuffle order
    do i=1, length*length
      rind = (rands(i) * length * length)
      if( rind .eq. 0) then
        rind = 1
      endif
      tmp = shuffled(i)
      shuffled(i) = shuffled(rind)
      shuffled(rind) = tmp
    end do
    trans = 0

    !> send back to pes
    do i=0, npes-1
      trans = RESHAPE(shuffled(SIZE(trans)*i + 1:SIZE(trans)*(i+1)), &
                         (/SIZE(trans,1), SIZE(trans,2) /) )
      if(i.ne.root) then
        call mpp_send(trans, SIZE(trans), i)
      else
        dataR4 = trans
      endif
    end do
  else
    call mpp_send(dataR4, SIZE(dataR4), root)
    call mpp_recv(trans, SIZE(dataR4), root)
    dataR4 = trans
  endif
  deallocate(trans, shuffled)
end subroutine shuffleDataR4

!> aggregates 64-bit real data on root and randomizes ordering, then sends partitions back to pes
subroutine shuffleDataR8(dataR8)
  real(r8_kind), intent(INOUT) :: dataR8(:,:)
  real(r8_kind), allocatable :: trans(:,:), shuffled(:), tmp
  integer :: rind

  allocate(trans(SIZE(dataR8,1), SIZE(dataR8,2)))
  allocate(shuffled(1:length*length))

  if( pe.eq.root) then
    !> get array partitions and aggregate into 1d
    shuffled(1:SIZE(dataR8)) = RESHAPE(dataR8, (/SIZE(dataR8)/))
    do i=1, npes-1
      call mpp_recv(trans, SIZE(dataR8) , i)
      shuffled( SIZE(trans)*i+1 : SIZE(trans)*(i+1)) = RESHAPE(trans, (/SIZE(trans)/))
    end do

    !> shuffle order
    do i=1, length*length
      rind = (rands(i) * length * length)
      if( rind .eq. 0) then
        rind = 1
      endif
      tmp = shuffled(i)
      shuffled(i) = shuffled(rind)
      shuffled(rind) = tmp
    end do
    trans = 0

    !> send back to pes
    do i=0, npes-1
      trans = RESHAPE(shuffled(SIZE(trans)*i + 1:SIZE(trans)*(i+1)), &
                         (/SIZE(trans,1), SIZE(trans,2) /) )
      if(i.ne.root) then
        call mpp_send(trans, SIZE(trans), i)
      else
        dataR8 = trans
      endif
    end do
  else
    call mpp_send(dataR8, SIZE(dataR8), root)
    call mpp_recv(trans, SIZE(dataR8), root)
    dataR8 = trans
  endif
  deallocate(trans, shuffled)
end subroutine shuffleDataR8

end program test_global_arrays
