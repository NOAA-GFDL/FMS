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
  use mpp_domains_mod, only: mpp_domains_init, mpp_define_domains, domain2d
  use mpp_domains_mod, only: mpp_define_layout, mpp_domains_set_stack_size
  use mpp_domains_mod, only: mpp_get_global_domain, mpp_global_max
  use mpp_domains_mod, only: mpp_global_min, mpp_get_data_domain,mpp_get_compute_domain
  use mpp_domains_mod, only: mpp_domains_exit, mpp_update_domains
  use mpp_domains_mod, only: mpp_get_domain_shift, mpp_global_sum
  use mpp_domains_mod, only: CYCLIC_GLOBAL_DOMAIN, NORTH, EAST, CENTER, CORNER, BITWISE_EXACT_SUM
  use mpp_mod,         only: MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED, mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use fms_mod,         only: check_nml_error, input_nml_file

  implicit none

  integer, parameter            :: length=64
  integer                       :: id, pe, npes, root, i, j, icount, jcount, io
  integer(i4_kind)              :: maxI4, minI4, ierr, sumI4, sumI4_5d, sumI4_shuf
  integer(i8_kind)              :: maxI8, minI8, sumI8, sumI8_5d, sumI8_shuf
  integer(i4_kind), allocatable :: dataI4(:,:), dataI4_shuf(:,:), recv_data_i4(:,:)
  integer(i8_kind), allocatable :: dataI8(:,:), dataI8_shuf(:,:), recv_data_i8(:,:)
  real(r4_kind), allocatable    :: dataR4(:,:), dataR4_shuf(:,:), recv_data_r4(:,:)
  real(r8_kind), allocatable    :: dataR8(:,:), dataR8_shuf(:,:), recv_data_r8(:,:)
  real, allocatable             :: rands(:)
  type(domain2D)                :: domain
  real(r8_kind)                 :: rcoef, maxR8, minR8, sumR8, sumR8_shuf
  real(r4_kind)                 :: maxR4, minR4, sumR4, sumR4_shuf
  integer                       :: isc, iec, jsc, jec
  integer                       :: isd, ied, jsd, jed
  character(len=32)             :: strTmp1, strTmp2
  integer(i4_kind), parameter   :: randmaxI4 = 2048
  integer(i8_kind), parameter   :: randmaxI8 = 4096
  real(r8_kind), parameter      :: tol4 = 1e-4, tol8 = 1e-6!> tolerance for real comparisons

  ! namelist variables - just logicals to enable individual tests
  ! simple just does normal max/min + sums across a domain
  ! full does max/min+sums with halos and symmetry
  logical :: test_simple= .false. , test_full = .false.
  namelist / test_global_arrays_nml / test_simple, test_full

  call mpp_init()

  call mpp_domains_init()
  !call mpp_set_stack_size(3145746)
  call mpp_domains_set_stack_size(4000000)

  read(input_nml_file, nml=test_global_arrays_nml, iostat=io)
  ierr = check_nml_error(io, 'test_global_arrays_nml')
  pe = mpp_pe()
  npes = mpp_npes()
  call mpp_set_root_pe(0)
  root = mpp_root_pe()
  if( test_simple) then
    call test_mpp_global_simple()
    deallocate(dataI4, dataI8, dataR4, dataR8, rands)
    deallocate(dataR4_shuf, dataR8_shuf,dataI4_shuf, dataI8_shuf)
  else if(test_full) then
    call test_global_reduce( 'Simple')
    call test_global_reduce( 'Simple symmetry center')
    call test_global_reduce( 'Simple symmetry corner')
    call test_global_reduce( 'Simple symmetry east')
    call test_global_reduce( 'Simple symmetry north')
    call test_global_reduce( 'Cyclic symmetry center')
    call test_global_reduce( 'Cyclic symmetry corner')
    call test_global_reduce( 'Cyclic symmetry east')
    call test_global_reduce( 'Cyclic symmetry north')
  else
    call mpp_error(FATAL, "test_global_arrays: either test_sum or test_max_min must be true in input.nml")
  endif
  call mpp_sync()

  call mpp_domains_exit()
  call MPI_FINALIZE(ierr)

  contains

subroutine test_mpp_global_simple()

  !> define domains and allocate
  call mpp_define_domains( (/1,length,1,length/), (/1,8/), domain, xhalo=0)
  call mpp_get_compute_domain(domain, jsc, jec, isc, iec)
  call mpp_get_data_domain(domain, jsd, jed, isd, ied)
  allocate(dataI4(jsd:jed, isd:ied),dataI8(jsd:jed, isd:ied), rands(length*length))
  allocate(dataR4(jsd:jed, isd:ied), dataR8(jsd:jed, isd:ied))
  allocate(dataR4_shuf(jsd:jed, isd:ied), dataR8_shuf(jsd:jed, isd:ied))
  allocate(dataI4_shuf(jsd:jed, isd:ied), dataI8_shuf(jsd:jed, isd:ied))
  allocate(recv_data_r4(jsd:jed, isd:ied), recv_data_r8(jsd:jed, isd:ied))
  allocate(recv_data_i4(jsd:jed, isd:ied), recv_data_i8(jsd:jed, isd:ied))

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
  call mpp_sync()

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

  !> moves the data into different pe's and checks the sum still matches
  dataR4_shuf = dataR4 ; dataR8_shuf = dataR8
  dataI4_shuf = dataI4 ; dataI8_shuf = dataI8
  !! swap data with neighboring pe
  if(modulo(pe, 2) .eq. 0) then
    print *, pe, pe+1, SUM(dataR8_shuf)
    call mpp_send(dataR4_shuf, SIZE(dataR4_shuf), pe+1)
    call mpp_recv(recv_data_r4, SIZE(dataR4_shuf), pe+1)
    call mpp_sync()
    call mpp_send(dataR8_shuf, SIZE(dataR8_shuf), pe+1)
    call mpp_recv(recv_data_r8, SIZE(dataR8_shuf), pe+1)
    call mpp_sync()
    call mpp_send(dataI4_shuf, SIZE(dataI4_shuf), pe+1)
    call mpp_recv(recv_data_I4, SIZE(dataI4_shuf), pe+1)
    call mpp_sync()
    call mpp_send(dataI8_shuf, SIZE(dataI8_shuf), pe+1)
    call mpp_recv(recv_data_I8, SIZE(dataI8_shuf), pe+1)
  else
    print *, pe, pe-1, SUM(dataR8_shuf)
    call mpp_recv(recv_data_r4, SIZE(dataR4_shuf), pe-1)
    call mpp_send(dataR4_shuf, SIZE(dataR4_shuf), pe-1)
    call mpp_sync()
    call mpp_recv(recv_data_r8, SIZE(dataR8_shuf), pe-1)
    call mpp_send(dataR8_shuf, SIZE(dataR8_shuf), pe-1)
    call mpp_sync()
    call mpp_send(dataI4_shuf, SIZE(dataI4_shuf), pe-1)
    call mpp_recv(recv_data_I4, SIZE(dataI4_shuf), pe-1)
    call mpp_sync()
    call mpp_send(dataI8_shuf, SIZE(dataI8_shuf), pe-1)
    call mpp_recv(recv_data_I8, SIZE(dataI8_shuf), pe-1)
  endif
  call mpp_sync()
  dataR4_shuf = recv_data_r4
  dataR8_shuf = recv_data_r8

  call mpp_error(NOTE, "----------Testing 32-bit real mpp_global_sum with reordering----------")
  call mpp_update_domains(dataR4_shuf, domain)
  sumR4_shuf = mpp_global_sum(domain, dataR4_shuf)
  ! check that shuffled array results are approximately the same as the original array
  if(abs(sumR4-sumR4_shuf) .gt. 1E-4 ) then
    strTmp1 = ""; strTmp2=""
    write(strTmp1,*) sumR4_shuf
    write(strTmp2,*) sumR4
    call mpp_error(FATAL,"test_global_arrays: invalid 32-bit real answer after reordering"// &
                   NEW_LINE('a')//"Sum: "// strTmp1// " ne "//strTmp2)
  endif

  call mpp_sync()
  call mpp_error(NOTE, "----------Testing 64-bit real mpp_global_sum with reordering----------")
  call mpp_update_domains(dataR8_shuf, domain)
  sumR8_shuf = mpp_global_sum(domain, dataR8_shuf)
  ! check that shuffled array results are approximately the same as the original array
  if(abs(sumR8-sumR8_shuf) .gt. 1E-7) then
    strTmp1 = ""; strTmp2=""
    write(strTmp1,*) sumR8_shuf
    write(strTmp2,*) sumR8
    call mpp_error(FATAL,"test_global_arrays: invalid 64-bit real answer after reordering"// &
                   NEW_LINE('a')//"Sum: "// strTmp1// " ne "//strTmp2)
  endif

  call mpp_error(NOTE, "----------Testing 32-bit integer mpp_global_sum with reordering----------")
  call mpp_update_domains(dataI4_shuf, domain)
  sumI4_shuf = mpp_global_sum(domain, dataI4_shuf)

  ! check that shuffled array results are approximately the same as the original array
  if(sumI4 .ne. sumI4_shuf) then
    strTmp1 = ""; strTmp2=""
    write(strTmp1,*) sumI4_shuf
    write(strTmp2,*) sumI4
    call mpp_error(FATAL,"test_global_arrays: invalid 32-bit integer answer after reordering"// &
                   NEW_LINE('a')//"Sum: "// strTmp1// " ne "//strTmp2)
  endif

  call mpp_error(NOTE, "----------Testing 64-bit integer mpp_global_sum with reordering----------")
  call mpp_update_domains(dataI8_shuf, domain)
  sumI8_shuf = mpp_global_sum(domain, dataI8_shuf)

  ! check that shuffled array results are approximately the same as the original array
  if(sumI8 .ne. sumI8_shuf) then
    strTmp1 = ""; strTmp2=""
    write(strTmp1,*) sumI8_shuf
    write(strTmp2,*) sumI8
    call mpp_error(FATAL,"test_global_arrays: invalid 64-bit integer answer after reordering"// &
                   NEW_LINE('a')//"Sum: "// strTmp1// " ne "//strTmp2)
  endif
end subroutine test_mpp_global_simple

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
    checkResultReal4 = checkResultReal4 .and. (abs((res(1)-tres(1))/res(1)) .lt. tol4) .and. &
                       (abs((res(2)-tres(2))/res(2)) .lt. tol4)
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
    checkResultReal8 = checkResultReal8 .and. (abs((res(1)-tres(1))/res(1)) .lt. tol8) .and. &
                       (abs((res(2)-tres(2))/res(2)) .lt. tol8)
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

  allocate(recv(2))
  ! root receives and sums local sums from each pe
  if(pe .eq. root) then
    nsum = SUM(dataR4)
    do i=1, npes - 1
      call mpp_recv(recv, 2, i)
      nsum = nsum + recv(1)
      ! also check for matching global sum
      if( abs((recv(2)-gsum)/gsum) .gt. tol4) then
        checkSumReal4 = .false.
        deallocate(recv)
        return
      endif
    end do
    checkSumReal4 = (abs((nsum-gsum)/gsum) .lt. tol4)
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

  allocate(recv(2))
  ! root receives and sums local sums from each pe
  if(pe .eq. root) then
    nsum = SUM(dataR8)
    do i=1, npes - 1
      call mpp_recv(recv, 2, i)
      nsum = nsum + recv(1)
      ! also check for matching global sum
      if( abs((recv(2)-gsum)/gsum) .gt. tol8) then
        checkSumReal8 = .false.
        deallocate(recv)
        return
      endif
    end do
    checkSumReal8 = (abs((nsum-gsum)/gsum) .lt. tol8)
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

   !--- test mpp_global_sum, mpp_global_min and mpp_global_max
  subroutine test_global_reduce (type)
    character(len=*), intent(in) :: type
    real    :: lsum, gsum, lmax, gmax, lmin, gmin
    integer :: ni, nj, ishift, jshift, position, k
    integer :: is, ie, js, je !, isd, ied, jsd, jed
    integer :: nx=128, ny=128, nz=40, stackmax=4000000
    integer :: layout(2)
    integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2
    real, allocatable, dimension(:,:,:) :: global1, x
    real, allocatable, dimension(:,:)   :: global2D
    !--- set up domain
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Simple' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                    shalo=shalo, nhalo=nhalo, name=type )
    case( 'Simple symmetry center', 'Simple symmetry corner', 'Simple symmetry east', 'Simple symmetry north' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                    shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case( 'Cyclic symmetry center', 'Cyclic symmetry corner', 'Cyclic symmetry east', 'Cyclic symmetry north' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo,&
                                    name=type, symmetry = .true., xflags=CYCLIC_GLOBAL_DOMAIN, &
                                            &  yflags=CYCLIC_GLOBAL_DOMAIN )
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !--- determine if an extra point is needed
    ishift = 0; jshift = 0; position = CENTER
    select case(type)
    case ('Simple symmetry corner', 'Cyclic symmetry corner')
       ishift = 1; jshift = 1; position = CORNER
    case ('Simple symmetry east', 'Cyclic symmetry east' )
       ishift = 1; jshift = 0; position = EAST
    case ('Simple symmetry north', 'Cyclic symmetry north')
       ishift = 0; jshift = 1; position = NORTH
    end select

    ie  = ie+ishift;  je  = je+jshift
    ied = ied+ishift; jed = jed+jshift
    ni  = nx+ishift;  nj  = ny+jshift
    allocate(global1(1-whalo:ni+ehalo, 1-shalo:nj+nhalo, nz))
    global1 = 0.0
    do k = 1,nz
       do j = 1,nj
          do i = 1,ni
             global1(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    enddo

    !--- NOTE: even though the domain is cyclic, no need to apply cyclic condition on the global data

    allocate( x (isd:ied,jsd:jed,nz) )
    allocate( global2D(ni,nj))

    x(:,:,:) = global1(isd:ied,jsd:jed,:)
    do j = 1, nj
       do i = 1, ni
          global2D(i,j) = sum(global1(i,j,:))
       enddo
    enddo
    !test mpp_global_sum

    if(type(1:6) == 'Simple') then
       gsum = sum( global2D(1:ni,1:nj) )
    else
       gsum = sum( global2D(1:nx, 1:ny) )
    endif
    id = mpp_clock_id( type//' sum', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lsum = mpp_global_sum( domain, x, position = position  )
    call mpp_clock_end  (id)
    if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', type//' Fast sum=', lsum, gsum

    !test exact mpp_global_sum
    id = mpp_clock_id( type//' exact sum', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lsum = mpp_global_sum( domain, x, BITWISE_EXACT_SUM, position = position )
    call mpp_clock_end  (id)
    !--- The following check will fail on altix in normal mode, but it is ok
    !--- in debugging mode. It is ok on irix.
    call compare_data_scalar(lsum, gsum, FATAL, type//' mpp_global_exact_sum')

    !test mpp_global_min
    gmin = minval(global1(1:ni, 1:nj, :))
    id = mpp_clock_id( type//' min', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lmin = mpp_global_min( domain, x, position = position )
    call mpp_clock_end  (id)
    call compare_data_scalar(lmin, gmin, FATAL, type//' mpp_global_min')

    !test mpp_global_max
    gmax = maxval(global1(1:ni, 1:nj, :))
    id = mpp_clock_id( type//' max', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lmax = mpp_global_max( domain, x, position = position )
    call mpp_clock_end  (id)
    call compare_data_scalar(lmax, gmax, FATAL, type//' mpp_global_max' )

    deallocate(global1, x)

  end subroutine test_global_reduce

  subroutine compare_data_scalar( a, b, action, string )
    real,             intent(in) :: a, b
    integer,          intent(in) :: action
    character(len=*), intent(in) :: string
    if( a .EQ. b)then
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': data comparison are OK.' )
    else
        call mpp_error( action, trim(string)//': data comparison are not OK.' )
    end if

  end subroutine compare_data_scalar

end program test_global_arrays
