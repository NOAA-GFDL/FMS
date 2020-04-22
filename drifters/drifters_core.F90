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
!
! nf95 -r8 -g -I ~/regression/ia64/23-Jun-2005/CM2.1U_Control-1990_E1.k32pe/include/ -D_TEST_DRIFTERS -D_F95 quicksort.F90 drifters_core.F90


module drifters_core_mod
#include <fms_platform.h>
  implicit none
  private

  public :: drifters_core_type, drifters_core_new, drifters_core_del, drifters_core_set_ids
  public :: drifters_core_remove_and_add, drifters_core_set_positions, assignment(=)
!#ifdef _TEST_DRIFTERS_CORE
  public :: drifters_core_print,  drifters_core_resize
!#endif

  ! Globals
  integer, parameter, private   :: MAX_STR_LEN = 128
! Include variable "version" to be written to log file.
#include<file_version.h>

  type drifters_core_type
     ! Be sure to update drifters_core_new, drifters_core_del and drifters_core_copy_new
     ! when adding members
     integer*8 :: it   ! time index
     integer :: nd     ! number of dimensions
     integer :: np     ! number of particles (drifters)
     integer :: npdim  ! max number of particles (drifters)
     integer, allocatable :: ids(:)_NULL  ! particle id number
     real   , allocatable :: positions(:,:)  
  end type drifters_core_type

  interface assignment(=)
     module procedure drifters_core_copy_new
  end interface

contains

!###############################################################################
  subroutine drifters_core_new(self, nd, npdim, ermesg)
    type(drifters_core_type)        :: self
    integer, intent(in)       :: nd
    integer, intent(in)       :: npdim
    character(*), intent(out) :: ermesg
    integer ier, iflag, i
    ermesg = ''
    ier    = 0

    call drifters_core_del(self, ermesg)

    allocate(self%positions(nd, npdim), stat=iflag)
    if(iflag/=0) ier = ier + 1
    self%positions   = 0.

    allocate(self%ids(npdim), stat=iflag)
    if(iflag/=0) ier = ier + 1
    self%ids         = (/(i, i=1,npdim)/)

    self%nd    = nd
    self%npdim = npdim

    if(ier/=0) ermesg = 'drifters::ERROR in drifters_core_new'
  end subroutine drifters_core_new

 !###############################################################################
 subroutine drifters_core_del(self, ermesg)
    type(drifters_core_type)        :: self
    character(*), intent(out) :: ermesg
    integer ier, iflag
    ermesg = ''
    ier    = 0
    self%it  = 0
    self%nd  = 0
    self%np  = 0
    iflag = 0
    if(allocated(self%positions)) deallocate(self%positions, stat=iflag)
    if(iflag/=0) ier = ier + 1
    if(allocated(self%ids)) deallocate(self%ids, stat=iflag)
    if(iflag/=0) ier = ier + 1

    if(ier/=0) ermesg = 'drifters::ERROR in drifters_core_del'
  end subroutine drifters_core_del

 !###############################################################################
 subroutine drifters_core_copy_new(new_instance, old_instance)

    type(drifters_core_type), intent(inout)   :: new_instance
    type(drifters_core_type), intent(in)      :: old_instance

    character(len=MAX_STR_LEN) :: ermesg

    ermesg = ''
    call drifters_core_del(new_instance, ermesg)
    if(ermesg/='') return
    ! this should provide the right behavior for both pointers and allocatables
    new_instance%it         = old_instance%it
    new_instance%nd         = old_instance%nd
    new_instance%np         = old_instance%np
    new_instance%npdim      = old_instance%npdim
    allocate(new_instance%ids( size(old_instance%ids) ))
    new_instance%ids        = old_instance%ids
    allocate(new_instance%positions( size(old_instance%positions,1), &
         &                           size(old_instance%positions,2) ))
    new_instance%positions  = old_instance%positions

 end subroutine drifters_core_copy_new
 !###############################################################################
  subroutine drifters_core_resize(self, npdim, ermesg)
    type(drifters_core_type)        :: self
    integer, intent(in)        :: npdim ! new max value
    character(*), intent(out) :: ermesg
    integer ier, iflag, i

    real   , allocatable :: positions(:,:)
    integer, allocatable :: ids(:)

    ermesg = ''
    ier    = 0
    if(npdim <= self%npdim) return

    ! temps
    allocate(positions(self%nd, self%np), stat=iflag)
    allocate(               ids(self%np), stat=iflag)

    positions    = self%positions(:, 1:self%np)
    ids          = self%ids(1:self%np)

    deallocate(self%positions, stat=iflag)
    deallocate(self%ids      , stat=iflag)

    allocate(self%positions(self%nd, npdim), stat=iflag)
    allocate(self%ids(npdim), stat=iflag)
    self%positions = 0.0
    ! default id numbers
    self%ids       = (/ (i, i=1,npdim) /)
    self%positions(:, 1:self%np) = positions
    self%npdim = npdim

    if(ier/=0) ermesg = 'drifters::ERROR in drifters_core_resize'
  end subroutine drifters_core_resize

!###############################################################################
  subroutine drifters_core_set_positions(self, positions, ermesg)
    type(drifters_core_type)        :: self
    real, intent(in)           :: positions(:,:)
    character(*), intent(out)  :: ermesg
    integer ier !, iflag
    ermesg = ''
    ier = 0
    self%np = min(self%npdim, size(positions, 2))
    self%positions(:,1:self%np) = positions(:,1:self%np)
    self%it                = self%it + 1
    if(ier/=0) ermesg = 'drifters::ERROR in drifters_core_set_positions'
  end subroutine drifters_core_set_positions

!###############################################################################
  subroutine drifters_core_set_ids(self1, ids1, ermesg1)
    type(drifters_core_type)        :: self1
    integer, intent(in)        :: ids1(:)
    character(*), intent(out)  :: ermesg1
    integer ier, np !, iflag
    ermesg1 = ''
    ier = 0
    np = min(self1%npdim, size(ids1))
    self1%ids(1:np) = ids1(1:np)
    if(ier/=0) ermesg1 = 'drifters::ERROR in drifters_core_set_ids'
  end subroutine drifters_core_set_ids

!###############################################################################
subroutine drifters_core_remove_and_add(self, indices_to_remove_in, &
     & ids_to_add, positions_to_add, &
     & ermesg)
    type(drifters_core_type)        :: self
    integer, intent(in   )     :: indices_to_remove_in(:)
    integer, intent(in   )     :: ids_to_add(:)
    real   , intent(in   )     :: positions_to_add(:,:)
    character(*), intent(out)  :: ermesg
    integer ier, np_add, np_remove, i, j, n_diff !, iflag
    integer indices_to_remove(size(indices_to_remove_in))
    external qksrt_quicksort
    ermesg = ''
    ier = 0

    ! copy, required so we can have indices_to_remove_in intent(in)
    indices_to_remove = indices_to_remove_in
    np_remove = size(indices_to_remove)
    np_add    = size(ids_to_add, 1)
    n_diff = np_add - np_remove

    ! cannot remove more than there are elements
    if(self%np + n_diff < 0) then
       ermesg = 'drifters::ERROR attempting to remove more elements than there are elements in drifters_core_remove_and_add'
       return
    endif

    ! check for overflow, and resize if necessary
    if(self%np + n_diff > self%npdim)  &
         & call drifters_core_resize(self, int(1.2*(self%np + n_diff))+1, ermesg)

    do i = 1, min(np_add, np_remove)
       j = indices_to_remove(i)
       self%ids(j)            = ids_to_add(i)
       self%positions(:,j)    = positions_to_add(:,i)
    enddo

    if(n_diff > 0) then
       ! all the particles to remove were removed and replaced. Just need to append
       ! remaining particles to end of list
       self%ids(         self%np+1:self%np+n_diff)   = ids_to_add(        np_remove+1:np_add)
       self%positions(:, self%np+1:self%np+n_diff)   = positions_to_add(:,np_remove+1:np_add)

       self%np = self%np + n_diff

    else if(n_diff < 0) then
       ! all the particles were added by filling in holes left by particles that
       ! were previously removed. Now remove remaining particles, starting from the end,
       ! by replacing the missing particle with a copy from the end.

       ! sort remaining indices in ascending order
       call qksrt_quicksort(size(indices_to_remove), indices_to_remove, np_add+1, np_remove)

       do i = np_remove, np_add+1, -1
          if(self%np <= 0) exit
          j = indices_to_remove(i)
          self%ids      (  j)    = self%ids      (  self%np)
          self%positions(:,j)    = self%positions(:,self%np)
          self%np = self%np - 1
       enddo
    endif

    if(ier/=0) ermesg = 'drifters::ERROR in drifters_core_remove_and_add'
  end subroutine drifters_core_remove_and_add

!###############################################################################
  subroutine drifters_core_print(self1, ermesg1)
    type(drifters_core_type)        :: self1
    character(*), intent(out) :: ermesg1
    integer j
    ermesg1 = ''

    print '(a,i10,a,i6,a,i6,a,i4,a,i4,a,i4)','it=',self1%it,  &
         & ' np=', self1%np, ' npdim=', self1%npdim

    print *,'ids and positions:'
    do j = 1, self1%np
       print *,self1%ids(j), self1%positions(:,j)
    enddo

  end subroutine drifters_core_print


end module drifters_core_mod
!###############################################################################
!###############################################################################
