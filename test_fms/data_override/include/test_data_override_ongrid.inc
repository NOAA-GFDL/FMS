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

subroutine COMPARE_DATA_ (Domain_in, actual_result, expected_result, p)
  integer, parameter                      :: lkind = TEST_FMS_KIND_ !< Real precision of the test
  type(domain2d), intent(in)              :: Domain_in             !< Domain with mask table
  real(lkind), dimension(:,:), intent(in) :: actual_result         !< Result from data_override
  real(lkind), intent(in)                 :: expected_result       !< Expected result from data_override
  integer, intent(in)                     :: p                     !< Dimension order permutation index
  integer                                 :: xsizec, ysizec        !< Size of the compute domain
  integer                                 :: xsized, ysized        !< Size of the data domain
  integer                                 :: n(2)                  !< Size of acual_result
  integer                                 :: nh(2)                 !< Size of the halos
  integer                                 :: i, j                  !< Helper indices

  !< Data is only expected to be overridden for the compute domain -not at the halos.
  call mpp_get_compute_domain(Domain_in, xsize=xsizec, ysize=ysizec)
  call mpp_get_data_domain(Domain_in, xsize=xsized, ysize=ysized)

  !< Note that actual_result has indices at (1:n(1),1:n(2)) not (is:ie,js:je)
  n(1) = size(actual_result, 1)
  n(2) = size(actual_result, 2)

  nh(1) = (xsized-xsizec)/2
  nh(2) = (ysized-ysizec)/2
  call permute_arr(nh, p)

  do i = 1, n(1)
    do j = 1, n(2)
      if (i <= nh(1) .or. i > (n(1)-nh(1)) .or. j <= nh(2) .or. j > (n(2)-nh(2))) then
        !< This is the result at the halos it should 999.
        if (actual_result(i,j) .ne. 999._lkind) then
          print *, "for i=", i, " and j=", j, " result=", actual_result(i,j)
          call mpp_error(FATAL, "test_data_override_ongrid: Data was overridden in the halos!!")
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
subroutine ONGRID_TEST_ (p)
  integer, intent(in) :: p
  integer, parameter                       :: lkind = TEST_FMS_KIND_ !< Real precision of the test
  real(lkind)                              :: expected_result        !< Expected result from data_override
  type(time_type)                          :: Time                   !< Time
  real(lkind), allocatable, dimension(:,:) :: runoff                 !< Data to be written
  type(permutable_indices(2))              :: ind                    !< Dimension boundaries

  ind%lb = [is, js]
  ind%ub = [ie, je]
  call ind%permute(p)

  allocate(runoff(ind%lb(1):ind%ub(1), ind%lb(2):ind%ub(2)))

  runoff = 999._lkind
  !< Run it when time=3
  Time = set_date(1,1,4,0,0,0)
  call data_override('OCN','runoff',runoff, Time)
  !< Because you are getting the data when time=3, and this is an "ongrid" case, the expected result is just
  !! equal to the data at time=3, which is 3.
  expected_result = 3._lkind
  call COMPARE_DATA_ (Domain, runoff, expected_result, p)

  !< Run it when time=4
  runoff = 999._lkind
  Time = set_date(1,1,5,0,0,0)
  call data_override('OCN','runoff',runoff, Time)
  !< You are getting the data when time=4, the data at time=3 is 3. and at time=5 is 4., so the expected result
  !! is the average of the 2 (because this is is an "ongrid" case and there is no horizontal interpolation).
  expected_result = (3._lkind + 4._lkind) / 2._lkind
  call COMPARE_DATA_ (Domain, runoff, expected_result, p)
end subroutine ONGRID_TEST_

!> @brief Tests bilinear data_override with and increasing and decreasing grid case
!! and comares the output betweeen the cases to ensure it is correct
subroutine BILINEAR_TEST_ (p)
  integer, intent(in) :: p
  integer, parameter                       :: lkind = TEST_FMS_KIND_ !< Real precision of the test
  type(time_type)                          :: Time                   !< Time
  real(lkind), allocatable, dimension(:,:) :: runoff_decreasing      !< Data to be written
  real(lkind), allocatable, dimension(:,:) :: runoff_increasing      !< Data to be written
  type(permutable_indices(2))              :: ind                    !< Dimension boundaries

  integer :: i, j, k
  logical :: success

  ind%lb = [is, js]
  ind%ub = [ie, je]
  call ind%permute(p)

  allocate(runoff_decreasing(ind%lb(1):ind%ub(1), ind%lb(2):ind%ub(2)))
  allocate(runoff_increasing(ind%lb(1):ind%ub(1), ind%lb(2):ind%ub(2)))

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
end subroutine BILINEAR_TEST_

subroutine WEIGHT_FILE_TEST_ (p)
  integer, intent(in) :: p
  integer, parameter                       :: lkind = TEST_FMS_KIND_ !< Real precision of the test
  type(time_type)                          :: Time                   !< Time
  real(lkind), allocatable, dimension(:,:) :: runoff                 !< Data from normal override
  real(lkind), allocatable, dimension(:,:) :: runoff_weight          !< Data from weight file override
  real(lkind)                              :: threshold              !< Threshold for the difference in answers
  type(permutable_indices(2))              :: ind                    !< Dimension boundaries

  integer :: i, j, k
  logical :: success

  ind%lb = [is, js]
  ind%ub = [ie, je]
  call ind%permute(p)

  allocate(runoff(ind%lb(1):ind%ub(1), ind%lb(2):ind%ub(2)))
  allocate(runoff_weight(ind%lb(1):ind%ub(1), ind%lb(2):ind%ub(2)))

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
end subroutine WEIGHT_FILE_TEST_

subroutine SCALAR_TEST_
  integer, parameter :: lkind = TEST_FMS_KIND_ !< Real precision of the test
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
  if (co2 .ne. expected_result) call mpp_error(FATAL, "co2 was not overridden to the correct value!")

  !< Run it when time=4
  co2 = 999._lkind
  Time = set_date(1,1,5,0,0,0)
  call data_override('OCN','co2',co2, Time)
  !< You are getting the data when time=4, the data at time=3 is 3. and at time=5 is 4., so the expected result
  !! is the average of the 2 (because this is is an "ongrid" case and there is no horizontal interpolation).
  expected_result = (3._lkind + 4._lkind) / 2._lkind
  if (co2 .ne. expected_result) call mpp_error(FATAL, "co2 was not overridden to the correct value!")
end subroutine SCALAR_TEST_

subroutine ENSEMBLE_TEST_ (p)
  integer, intent(in) :: p
  integer, parameter                       :: lkind = TEST_FMS_KIND_ !< Real precision of the test
  real(lkind)                              :: expected_result        !< Expected result from data_override
  type(time_type)                          :: Time                   !< Time
  real(lkind), allocatable, dimension(:,:) :: runoff                 !< Data to be written
  type(permutable_indices(2))              :: ind                    !< Dimension boundaries
  integer                                  :: scale_fac              !< Scale factor to use when determining
                                                                     !! the expected answer
  logical :: sucessful !< .True. if the data_override was sucessful

  ind%lb = [is, js]
  ind%ub = [ie, je]
  call ind%permute(p)

  allocate(runoff(ind%lb(1):ind%ub(1), ind%lb(2):ind%ub(2)))

  scale_fac = ensemble_id
  if (test_case .eq. ensemble_same_yaml) scale_fac = 1

  runoff = 999._lkind
  !< Run it when time=3
  Time = set_date(1,1,4,0,0,0)
  call data_override('OCN','runoff',runoff, Time, override=sucessful)
  if (.not. sucessful) call mpp_error(FATAL, "The data was not overridden correctly")
  !< Because you are getting the data when time=3, and this is an "ongrid" case, the expected result is just
  !! equal to the data at time=3, which is 3+scale_fac.
  expected_result = 3._lkind + real(scale_fac,kind=lkind)
  call COMPARE_DATA_ (Domain, runoff, expected_result, p)

  !< Run it when time=4
  runoff = 999._lkind
  Time = set_date(1,1,5,0,0,0)
  call data_override('OCN','runoff',runoff, Time, override=sucessful)
  if (.not. sucessful) call mpp_error(FATAL, "The data was not overridden correctly")
  !< You are getting the data when time=4, the data at time=3 is 3+scale_fac. and at time=5 is 4+scale_fac.,
  !! so the expected result is the average of the 2 (because this is is an "ongrid" case and there
  !! is no horizontal interpolation).
  expected_result = (3._lkind + real(scale_fac,kind=lkind) + 4._lkind + real(scale_fac,kind=lkind)) / 2._lkind
  call COMPARE_DATA_ (Domain, runoff, expected_result, p)
end subroutine ENSEMBLE_TEST_

subroutine MULTI_FILE_TESTS_ (p)
  integer, intent(in) :: p
  integer, parameter                       :: lkind = TEST_FMS_KIND_ !< Real precision of the test
  real(lkind)                              :: expected_result        !< Expected result from data_override
  type(time_type)                          :: Time                   !< Time
  real(lkind), allocatable, dimension(:,:) :: runoff                 !< Data to be written
  type(permutable_indices(2))              :: ind                    !< Dimension boundaries
  logical :: sucessful !< .True. if the data_override was sucessful

  ind%lb = [is, js]
  ind%ub = [ie, je]
  call ind%permute(p)

  allocate(runoff(ind%lb(1):ind%ub(1), ind%lb(2):ind%ub(2)))

  !< Run it when time=3, this is going to use the previous file and the current file
  Time = set_date(1,1,4,0,0,0)
  runoff = 999._lkind
  call data_override('OCN','runoff',runoff, Time, override=sucessful)
  if (.not. sucessful) call mpp_error(FATAL, "The data was not overridden correctly")
  expected_result = 3._lkind
  call COMPARE_DATA_ (Domain, runoff, expected_result, p)

  !< Run it when time=4, this is going to use the current file
  Time = set_date(1,1,5,0,0,0)
  runoff = 999._lkind
  call data_override('OCN','runoff',runoff, Time, override=sucessful)
  if (.not. sucessful) call mpp_error(FATAL, "The data was not overridden correctly")
  expected_result = 4._lkind
  call COMPARE_DATA_ (Domain, runoff, expected_result, p)

  !< Run it when time=5, this is going to use the current file and the next file
  Time = set_date(1,1,6,0,0,0)
  runoff = 999._lkind
  call data_override('OCN','runoff',runoff, Time, override=sucessful)
  if (.not. sucessful) call mpp_error(FATAL, "The data was not overridden correctly")
  expected_result = 5._lkind
  call COMPARE_DATA_ (Domain, runoff, expected_result, p)
end subroutine MULTI_FILE_TESTS_
