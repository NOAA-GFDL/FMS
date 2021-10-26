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
    call test_reproducing_sum_r8()
  case(2)
    call test_efp_type()
  case default
    call mpp_error(FATAL, "invalid test_num given")
  end select

  call mpp_exit() 

  contains

  subroutine test_reproducing_sum_r8()
    real(r8_kind)        :: arr1(64,64), arr2(64,64), res, res_chk
    integer              :: err, i, j, n, npes
    type(mpp_efp_type)   :: res_efp
    integer, allocatable :: pelist(:)

    ! check normal functionality with same pe count
    arr1 = 0; arr2 = 0; res_chk = 0
    do i=1, 64
      do j=1, 64
        arr1(i,j) = i * 1d5 + j
        arr2(j,i) = i * 1d5 + j
        res_chk = res_chk + i * 1d5 + j
      end do
    end do
    res = mpp_reproducing_sum(arr1)
    if (res .ne. mpp_reproducing_sum(arr2) .or. res .ne. res_chk*mpp_npes()) call mpp_error(FATAL, &
                                        "test_reproducing_sum: reproducing sums differ")

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
    arr1(1,1) = 0.0 / 0.0 ! NaN
    res = mpp_reproducing_sum(arr1, overflow_check = .true., err = err)
    if ( res .ne. 0.0_r8_kind .or. err .ne. 4 ) call mpp_error(FATAL, "test_reproducing_sum: " &
                                                //"NaN error not detected")
    arr1(1,1) = 1.0 / 0.0 ! inf 
    res = mpp_reproducing_sum(arr1, overflow_check = .true., err = err)
    if ( res .ne. 0.0_r8_kind .or. err .ne. 2 ) call mpp_error(FATAL, "test_reproducing_sum: " &
                                                //"inf error not detected")
    arr1(1,1) = -1.0 / 0.0 ! -inf
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

  ! tests operations with mpp_efp_type
  subroutine test_efp_type()
    type(mpp_efp_type) :: efp1, efp2, efp_sum, efp_dif
    real :: rsum, rval1, rval2, rdif
    integer :: ex = 64

    rval1 = 2.**ex; rval2 = 2.**(ex-1)
    efp1 = mpp_real_to_efp(rval1)
    efp2 = mpp_real_to_efp(rval2)
    efp_sum = efp1 + efp2
    efp_dif = efp1 - efp2
    rsum = mpp_efp_to_real(efp_sum)
    rdif = mpp_efp_to_real(efp_dif)
    print *, rval1, rval2, rsum, rdif
    if( rsum .ne. 3*rval2 .or. rdif .ne. rval2) call mpp_error(FATAL, "test_efp_type: incorrect sum/difference") 

  end subroutine

end program
