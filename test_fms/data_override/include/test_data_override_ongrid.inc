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

subroutine COMPARE_DATA_ (Domain_in, actual_result, expected_result)
  integer, parameter                      :: lkind = DO_TEST_KIND_ !< Real precision of the test
  type(domain2d), intent(in)              :: Domain_in             !< Domain with mask table
  real(lkind), intent(in)                 :: expected_result       !< Expected result from data_override
  real(lkind), dimension(:,:), intent(in) :: actual_result         !< Result from data_override
  integer                                 :: xsizec, ysizec        !< Size of the compute domain
  integer                                 :: xsized, ysized        !< Size of the data domain
  integer                                 :: nx, ny                !< Size of acual_result
  integer                                 :: nhx, nhy              !< Size of the halos
  integer                                 :: i, j                  !< Helper indices

  !< Data is only expected to be overriden for the compute domain -not at the halos.
  call mpp_get_compute_domain(Domain_in, xsize=xsizec, ysize=ysizec)
  call mpp_get_data_domain(Domain_in, xsize=xsized, ysize=ysized)

  !< Note that actual_result has indices at (1:nx,1:ny) not (is:ie,js:je)
  nhx= (xsized-xsizec)/2
  nhy = (ysized-ysizec)/2
  nx = size(actual_result, 1)
  ny = size(actual_result, 2)

  do i = 1, nx
    do j = 1, ny
      if (i <= nhx .or. i > (nx-nhx) .or. j <= nhy .or. j > (ny-nhy)) then
        !< This is the result at the halos it should 999.
        if (actual_result(i,j) .ne. 999._lkind) then
          print *, "for i=", i, " and j=", j, " result=", actual_result(i,j)
          call mpp_error(FATAL, "test_data_override_ongrid: Data was overriden in the halos!!")
        endif
      else
        if (actual_result(i,j) .ne. expected_result) then
          print *, "for i=", i, " and j=", j, " result=", actual_result(i,j), " expected=", expected_result
          call mpp_error(FATAL, "test_data_override_ongrid: Result is different from expected answer!")
        endif
      endif
    enddo
  enddo
end subroutine COMPARE_DATA_

!> @brief Tests ongrid data overrides.
!! In the first case there is no time interpolation
!! In the second case there is time interpolation
subroutine ONGRID_TEST_
  integer, parameter                       :: lkind = DO_TEST_KIND_ !< Real precision of the test
  real(lkind)                              :: expected_result       !< Expected result from data_override
  type(time_type)                          :: Time                  !< Time
  real(lkind), allocatable, dimension(:,:) :: runoff                !< Data to be written

  allocate(runoff(is:ie,js:je))

  runoff = 999._lkind
  !< Run it when time=3
  Time = set_date(1,1,4,0,0,0)
  call data_override('OCN','runoff',runoff, Time)
  !< Because you are getting the data when time=3, and this is an "ongrid" case, the expected result is just
  !! equal to the data at time=3, which is 3.
  expected_result = 3._lkind
  call COMPARE_DATA_ (Domain, runoff, expected_result)

  !< Run it when time=4
  runoff = 999._lkind
  Time = set_date(1,1,5,0,0,0)
  call data_override('OCN','runoff',runoff, Time)
  !< You are getting the data when time=4, the data at time=3 is 3. and at time=5 is 4., so the expected result
  !! is the average of the 2 (because this is is an "ongrid" case and there is no horizontal interpolation).
  expected_result = (3._lkind + 4._lkind) / 2._lkind
  call COMPARE_DATA_ (Domain, runoff, expected_result)

  deallocate(runoff)
end subroutine ONGRID_TEST_

!> @brief Tests bilinear data_override with and increasing and decreasing grid case
!! and comares the output betweeen the cases to ensure it is correct
subroutine BILINEAR_TEST_
  integer, parameter                       :: lkind = DO_TEST_KIND_ !< Real precision of the test
  type(time_type)                          :: Time              !< Time
  real(lkind), allocatable, dimension(:,:) :: runoff_decreasing !< Data to be written
  real(lkind), allocatable, dimension(:,:) :: runoff_increasing !< Data to be written

  integer :: i, j, k
  logical :: success

  allocate(runoff_decreasing(is:ie,js:je))
  allocate(runoff_increasing(is:ie,js:je))

  runoff_decreasing = 999_lkind
  runoff_increasing = 999_lkind
  Time = set_date(1,1,4,0,0,0)
  call data_override('OCN','runoff_increasing',runoff_increasing, Time, override=success)
  if (.not. success) call mpp_error(FATAL, "Data override failed")
  call data_override('OCN','runoff_decreasing',runoff_decreasing, Time, override=success)
  if (.not. success) call mpp_error(FATAL, "Data override failed")

  do i = is, ie
    do j =  js, je
      if (abs(runoff_decreasing(i,j) - runoff_increasing(i,j)) .gt. 1) then
        call mpp_error(FATAL, "The data is not the same: "// &
        string(i)//","//string(j)//":"// &
        string(runoff_decreasing(i,j))//" vs "//string(runoff_increasing(i,j)))
      endif
    enddo
  enddo
  deallocate(runoff_decreasing, runoff_increasing)
end subroutine BILINEAR_TEST_

subroutine WEIGHT_FILE_TEST_
  integer, parameter                       :: lkind = DO_TEST_KIND_ !< Real precision of the test
  type(time_type)                          :: Time                  !< Time
  real(lkind), allocatable, dimension(:,:) :: runoff                !< Data from normal override
  real(lkind), allocatable, dimension(:,:) :: runoff_weight         !< Data from weight file override
  real(lkind)                              :: threshold             !< Threshold for the difference in answers

  integer :: i, j, k
  logical :: success

  allocate(runoff(is:ie,js:je))
  allocate(runoff_weight(is:ie,js:je))

  runoff = 999_lkind
  runoff_weight = 999_lkind
  Time = set_date(1,1,4,0,0,0)
  call data_override('OCN','runoff_obs',runoff, Time, override=success)
  if (.not. success) call mpp_error(FATAL, "Data override failed")
  call data_override('OCN','runoff_obs_weights',runoff_weight, Time, override=success)
  if (.not. success) call mpp_error(FATAL, "Data override failed")

  threshold = 1e-09
  if (lkind .eq. 4) then
    threshold = 1e-03
  endif

  do i = is, ie
    do j =  js, je
      if (abs(runoff(i,j) - runoff_weight(i,j)) .gt. threshold) then
        call mpp_error(FATAL, "The data is not the same: "// &
        string(i)//","//string(j)//":"// &
        string(runoff(i,j))//" vs "//string(runoff_weight(i,j)))
      endif
    enddo
  enddo
  deallocate(runoff, runoff_weight)
end subroutine WEIGHT_FILE_TEST_

subroutine SCALAR_TEST_
  integer, parameter :: lkind = DO_TEST_KIND_ !< Real precision of the test
  real(lkind)        :: expected_result       !< Expected result from data_override
  type(time_type)    :: Time                  !< Time
  real(lkind)        :: co2                   !< Data to be written

  co2 = 999._lkind
  !< Run it when time=3
  Time = set_date(1,1,4,0,0,0)
  call data_override('OCN','co2',co2, Time)
  !< Because you are getting the data when time=3, and this is an "ongrid" case, the expected result is just
  !! equal to the data at time=3, which is 3.
  expected_result = 3._lkind
  if (co2 .ne. expected_result) call mpp_error(FATAL, "co2 was not overriden to the correct value!")

  !< Run it when time=4
  co2 = 999._lkind
  Time = set_date(1,1,5,0,0,0)
  call data_override('OCN','co2',co2, Time)
  !< You are getting the data when time=4, the data at time=3 is 3. and at time=5 is 4., so the expected result
  !! is the average of the 2 (because this is is an "ongrid" case and there is no horizontal interpolation).
  expected_result = (3._lkind + 4._lkind) / 2._lkind
  if (co2 .ne. expected_result) call mpp_error(FATAL, "co2 was not overriden to the correct value!")
end subroutine SCALAR_TEST_

subroutine ENSEMBLE_TEST_
  integer, parameter                       :: lkind = DO_TEST_KIND_ !< Real precision of the test
  real(lkind)                              :: expected_result       !< Expected result from data_override
  type(time_type)                          :: Time                  !< Time
  real(lkind), allocatable, dimension(:,:) :: runoff                !< Data to be written
  integer                                  :: scale_fac             !< Scale factor to use when determining
                                                                    !! the expected answer
  logical :: sucessful !< .True. if the data_override was sucessful

  allocate(runoff(is:ie,js:je))

  scale_fac = ensemble_id
  if (test_case .eq. ensemble_same_yaml) scale_fac = 1

  runoff = 999._lkind
  !< Run it when time=3
  Time = set_date(1,1,4,0,0,0)
  call data_override('OCN','runoff',runoff, Time, override=sucessful)
  if (.not. sucessful) call mpp_error(FATAL, "The data was not overriden correctly")
  !< Because you are getting the data when time=3, and this is an "ongrid" case, the expected result is just
  !! equal to the data at time=3, which is 3+scale_fac.
  expected_result = 3._lkind + real(scale_fac,kind=lkind)
  call COMPARE_DATA_ (Domain, runoff, expected_result)

  !< Run it when time=4
  runoff = 999._lkind
  Time = set_date(1,1,5,0,0,0)
  call data_override('OCN','runoff',runoff, Time, override=sucessful)
  if (.not. sucessful) call mpp_error(FATAL, "The data was not overriden correctly")
  !< You are getting the data when time=4, the data at time=3 is 3+scale_fac. and at time=5 is 4+scale_fac.,
  !! so the expected result is the average of the 2 (because this is is an "ongrid" case and there
  !! is no horizontal interpolation).
  expected_result = (3._lkind + real(scale_fac,kind=lkind) + 4._lkind + real(scale_fac,kind=lkind)) / 2._lkind
  call COMPARE_DATA_ (Domain, runoff, expected_result)

  deallocate(runoff)
end subroutine ENSEMBLE_TEST_
