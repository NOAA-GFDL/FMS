! -*-f90-*-

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
!> @file
!> @brief Utility routines for parallelization, non-mpi version

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!         MISCELLANEOUS UTILITIES: mpp_error                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mpp_error_basic( errortype, errormsg )
  !a very basic error handler
  !uses ABORT and FLUSH calls, may need to use cpp to rename
  integer,                    intent(in) :: errortype
  character(len=*), intent(in), optional :: errormsg
  character(len=512)                     :: text
  logical                                :: opened
  integer                                :: istat, errunit, outunit

  if( .NOT.module_is_initialized )call ABORT()

  select case( errortype )
  case(NOTE)
     text = 'NOTE'         !just FYI
  case(WARNING)
     text = 'WARNING'      !probable error
  case(FATAL)
     text = 'FATAL'        !fatal error
  case default
     text = 'WARNING: non-existent errortype (must be NOTE|WARNING|FATAL)'
  end select

  if( npes.GT.1 )write( text,'(a,i5)' )trim(text)//' from PE', pe   !this is the mpp part
  if( PRESENT(errormsg) )text = trim(text)//': '//trim(errormsg)

  errunit = stderr()
  outunit = stdout()

  select case( errortype )
  case(NOTE)
     write( outunit,'(a)' )trim(text)
  case default
     write( errunit,'(/a/)' )trim(text)
     write( outunit,'(/a/)' )trim(text)
     if( errortype.EQ.FATAL .OR. warnings_are_fatal )then
        FLUSH(outunit)
        call ABORT() !automatically calls traceback on Cray systems
     end if
  end select

  error_state = errortype
  return
end subroutine mpp_error_basic

!#####################################################################
!> Makes a PE set out of a PE list. A PE list is an ordered list of PEs
!! a PE set is a triad (start,log2stride,size) for SHMEM, an a communicator for MPI
!! if stride is non-uniform or not a power of 2,
!! will return error (not required for MPI but enforced for uniformity)
function get_peset(pelist)
  integer                       :: get_peset
  integer, intent(in), optional :: pelist(:)

  if( .NOT.PRESENT(pelist) )then !set it to current_peset_num
     get_peset = current_peset_num; return
  end if

  get_peset = 0

  return

end function get_peset

!#######################################################################
!> synchronize PEs in list
subroutine mpp_sync( pelist, check )
  integer, intent(in), optional :: pelist(:)
  integer, intent(in), optional :: check

  return
end subroutine mpp_sync

!#######################################################################
!> This is to check if current PE's outstanding puts are complete
!! but we can't use shmem_fence because we are actually waiting for
!! a remote PE to complete its get
subroutine mpp_sync_self( pelist, check, request, msg_size, msg_type )
  integer, intent(in), optional :: pelist(:)
  integer, intent(in), optional :: check
  integer, intent(inout), optional :: request(:)
  integer, intent(in   ), optional :: msg_size(:)
  integer, intent(in   ), optional :: msg_type(:)


  return
end subroutine mpp_sync_self

function mpp_commID()
  integer :: mpp_commID

  call mpp_error(NOTE, "mpp_commID() is deprecated. Please use mpp_comm() instead.")
  mpp_commID = -1
end function mpp_commID

!> @brief Declare a pelist.
!!
!> This call is written specifically to accommodate a MPI restriction
!! that requires a parent communicator to create a child communicator, In
!! other words: a pelist cannot go off and declare a communicator, but
!! every PE in the parent, including those not in pelist(:), must get
!! together for the <TT>MPI_COMM_CREATE</TT> call. The parent is
!! typically <TT>MPI_COMM_WORLD</TT>, though it could also be a subset
!! that includes all PEs in <TT>pelist</TT>.
!!
!! This call implies synchronization across the PEs in the current
!! pelist, of which <TT>pelist</TT> is a subset.
subroutine mpp_declare_pelist(pelist, name, commID)
  integer,                    intent(in) :: pelist(:) !> pelist you are declaring and storing within FMS
  character(len=*), intent(in), optional :: name      !> unique name for an input pelist
  integer, intent(out), optional  :: commID           !> return of current MPI comm group communicator ID (deprecated)
  integer :: i

  if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DECLARE_PELIST: You must first call mpp_init.' )
  i = get_peset(pelist)
  write( peset(i)%name,'(a,i2.2)' ) 'PElist', i !default name
  if( PRESENT(name) )peset(i)%name = name
  if( PRESENT(commID) ) commID = -1
end subroutine mpp_declare_pelist

!this is created for use by mpp_define_domains within a pelist
!will be published but not publicized
subroutine mpp_get_current_pelist(pelist, name, commID)
  integer, intent(out) :: pelist(:)
  character(len=*), intent(out), optional :: name
  integer, intent(out), optional :: commID

  if( size(pelist(:)).NE.size(peset(current_peset_num)%list(:)) ) &
       call mpp_error( FATAL, 'MPP_GET_CURRENT_PELIST: size(pelist) is wrong.' )
  pelist(:) = peset(current_peset_num)%list(:)
  if( PRESENT(name) )name = peset(current_peset_num)%name
  if( PRESENT(commID) ) commID = -1
end subroutine mpp_get_current_pelist
