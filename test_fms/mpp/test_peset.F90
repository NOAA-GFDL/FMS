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
program test_peset
  !> @author Jessica Liptak
  !> @description Indirectly test the functionality of get_peset and expand_peset in mpp_util_sma.inc by passing
  !! processor IDs among PEs via mpp_send and mpp_recv, then calling mpp_sync_self.
  use mpp_mod, only : mpp_init, mpp_sync_self, mpp_exit, mpp_npes, &
                      mpp_error, FATAL, mpp_pe, mpp_init_test_peset_allocated, &
                      mpp_broadcast, mpp_root_pe, mpp_declare_pelist, input_nml_file
  implicit none

  integer :: npes !< The total number of PEs returned from mpp_npes
  integer, dimension(:), allocatable :: pelist !< List of pes that will communicate
  integer, dimension(:), allocatable :: numlist !< List of integers assigned to each pe
  integer :: i, ierr, sum_id
  !> Initialize mpp
  call mpp_init(test_level=mpp_init_test_peset_allocated)
  !> Get the total number of PEs
  npes = mpp_npes()
  !> create a pelist and a list of numbers to populate
  allocate(pelist(npes))
  allocate(numlist(npes))
  pelist=(/(i,i=0,npes-1)/)
  call mpp_declare_pelist(pelist)
  ! sum used to verify that all PEs received all numbers
  sum_id=sum(pelist)
  numlist(:)=0
  !> loop through number of pes and populate number list if the current pe is the root pe
  do i = 1, npes
    if (mpp_pe() .eq. mpp_root_pe()) &
      numlist(i) = i-1
  enddo
  !> root pe sends the numlist to all of the other pes
  call mpp_broadcast(numlist, npes, mpp_root_pe())
  !> check that calling routines are complete
  !! @note: mpp_sync_self calls get_peset; get pe_set calls expand_peset
  call mpp_sync_self(pelist)
  if (sum(numlist) .ne. sum_id) call mpp_error(FATAL,"numlist sum does not match pelist sum")
  if (allocated(numlist)) deallocate(numlist)
  if (allocated(pelist)) deallocate(pelist)
  !> Finalize mpp
  call MPI_FINALIZE(ierr)
end program test_peset
