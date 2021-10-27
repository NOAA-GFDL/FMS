program test_mpp_efp

  use fms

  implicit none

  integer :: test_num = 0, io_status
  namelist / test_mpp_efp_nml / test_num

  call mpp_init()

  read (input_nml_file, test_mpp_efp_nml, iostat=io_status)
  if (io_status > 0) call mpp_error(FATAL,'=>test_mpp_efp: Error reading input.nml')

  select case ( test_num )
  case(1)
    call test_reproducing_sum()
  case(2)
    call test_efp_type()
  case(3)
    call test_list_sum()
  case default
    call mpp_error(FATAL, "invalid test_num given")
  end select

  call mpp_exit()

  contains

  subroutine test_reproducing_sum()
    real(r8_kind)        :: arr1(64,64), arr2(64,64), arr3d(64,64,1), res, res_chk
    real(r4_kind)        :: arr_r4(64,64), res4
    integer              :: err, i, j, n, npes, zero
    type(mpp_efp_type)   :: res_efp
    integer, allocatable :: pelist(:)

    ! check normal functionality
    arr1 = 0; arr2 = 0; arr3d = 0; arr_r4 = 0
    res_chk = 0
    do i=1, 64
      do j=1, 64
        res_chk = res_chk + i * 1d5 + j
        arr1(i,j) = i * 1d5 + j
        arr2(j,i) = i * 1d5 + j
        arr3d(i, j, 1) = i * 1d5 + j
        arr_r4(i,j) = i * 1e5 + j
      end do
    end do
    res = mpp_reproducing_sum(arr1)
    if (res .ne. res_chk*mpp_npes()) call mpp_error(FATAL, "test_reproducing_sum: incorrect 2d sum")
    if (res .ne. mpp_reproducing_sum(arr2)) call mpp_error(FATAL, "test_reproducing_sum: 2d reproducing sums differ")
    if (res .ne. mpp_reproducing_sum(arr3d)) call mpp_error(FATAL, "test_reproducing_sum: 3d reproducing sums differ")
    res4 = mpp_reproducing_sum(arr_r4, efp_sum = res_efp) ! r4 needs to get checked through higher precision efp_type result
    if(res .ne. mpp_efp_to_real(res_efp)) call mpp_error(FATAL, "test_reproducing_sum: 2D r4 sum differs")

    ! test reproducibility across different pe counts
    npes = mpp_npes()
    do n=0, npes - 1
      allocate(pelist(n+1))
      pelist = (/( i, i = 0, n )/)
      call mpp_declare_pelist(pelist)
      if( any( pelist == mpp_pe() )) then
        call mpp_set_current_pelist(pelist)
        ! sum/npes should equal initial sum/all_npes
        if( mpp_reproducing_sum(arr1)/(n+1) .ne. res/npes) call mpp_error(FATAL, &
                               "test_reproducing_sum: reproducing sums differ between pe count")
      endif
      call mpp_set_current_pelist()
      deallocate(pelist)
    end do

    ! check edge cases
    zero = 0.0 ! avoid div by zero errors
    arr1(1,1) = 0.0 / zero ! NaN
    res = mpp_reproducing_sum(arr1, overflow_check = .true., err = err)
    if ( res .ne. 0.0_r8_kind .or. err .ne. 4 ) call mpp_error(FATAL, "test_reproducing_sum: " &
                                                //"NaN error not detected")
    arr1(1,1) = 1.0 / zero ! inf
    res = mpp_reproducing_sum(arr1, overflow_check = .true., err = err)
    if ( res .ne. 0.0_r8_kind .or. err .ne. 2 ) call mpp_error(FATAL, "test_reproducing_sum: " &
                                                //"inf error not detected")
    arr1(1,1) = -1.0 / zero ! -inf
    res = mpp_reproducing_sum(arr1, overflow_check = .true., err = err)
    if ( res .ne. 0.0_r8_kind .or. err .ne. 2 ) call mpp_error(FATAL, "test_reproducing_sum: " &
                                                //"-inf error not detected")
    arr1(1,1) = 1d47  ! single val overflow
    res = mpp_reproducing_sum(arr1, overflow_check = .true., err = err)
    if ( res .ne. 0.0_r8_kind .or. err .ne. 2 ) call mpp_error(FATAL, "test_reproducing_sum: " &
                                                //"single value overflow error not detected")
    arr1 = 1d45 ! sum overflow
    res = mpp_reproducing_sum(arr1, overflow_check = .true., err = err)
    if ( res .ne. 0.0_r8_kind .or. err .ne. 2 ) call mpp_error(FATAL, "test_reproducing_sum: " &
                                                //"sum value overflow error not detected")

  end subroutine

  ! tests operations overrides and conversions with mpp_efp_type
  subroutine test_efp_type()
    type(mpp_efp_type) :: efp1, efp2, efp_sum, efp_dif
    real(r8_kind) :: rsum, rval1, rval2, rdif
    integer :: ex = 64

    rval1 = 2.0**ex; rval2 = 2.0**(ex-1)
    efp1 = mpp_real_to_efp(rval1)
    efp2 = mpp_real_to_efp(rval2)
    efp_sum = efp1 + efp2
    efp_dif = efp1 - efp2
    rsum = mpp_efp_to_real(efp_sum)
    rdif = mpp_efp_to_real(efp_dif)
    if( rsum .ne. 3*rval2 .or. rdif .ne. rval2) call mpp_error(FATAL, "test_efp_type: "// &
                                         "incorrect sum/difference from operator overrides")
    if( rdif .ne. mpp_efp_real_diff(efp1, efp2))call mpp_error(FATAL, "test_efp_type: "// &
                                         "incorrect difference from mpp_efp_real_diff")

  end subroutine

  subroutine test_list_sum()
    type(mpp_efp_type) :: efps(64)
    integer :: i
    logical :: errors(64)
    real(r8_kind) :: rsum
    ! normal functionality
    do i=1, 64
      efps(i) = mpp_real_to_efp(2.0_r8_kind**i)
    end do
    call mpp_efp_list_sum_across_pes(efps, 64,errors)
    if(any(errors)) call mpp_error(FATAL, "test_list_sum: error from mpp_efp_list_sum_across_pes")
    do i=1, 64
      rsum = mpp_efp_to_real(efps(i))
      if(rsum .ne. (2.0_r8_kind**i) * mpp_npes()) call mpp_error(FATAL, "test_list_sum: "&
                              //"mpp_efp_list_sum_across_pes returned incorrect result")
    enddo
    ! catch overflow
    efps(64) = mpp_real_to_efp(1.e41)
    do i=1, 10
      call mpp_efp_list_sum_across_pes(efps, 64, errors)
    enddo
    if(.not. errors(64)) call mpp_error(FATAL, "test_list_sum: mpp_efp_list_sum_across_pes missed overflow error")

  end subroutine

end program
