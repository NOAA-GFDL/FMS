
! -*-f90-*-
! $Id: mpp_util_nocomm.inc,v 13.0 2006/03/28 21:42:00 fms Exp $

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
  character(len=256)                     :: text
  logical                                :: opened
  integer                                :: istat

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

  select case( errortype )
  case(NOTE)
     write( stdout(),'(a)' )trim(text)
  case default
     write( stderr(),'(/a/)' )trim(text)
     write( stdout(),'(/a/)' )trim(text)
     if( errortype.EQ.FATAL .OR. warnings_are_fatal )then
        call FLUSH(stdout(),istat)
#ifdef sgi_mipspro
        call TRACE_BACK_STACK_AND_PRINT()
#endif
        call ABORT() !automatically calls traceback on Cray systems
     end if
  end select

  error_state = errortype
  return
end subroutine mpp_error_basic

!#####################################################################
!overloads to mpp_error_basic
!support for error_mesg routine in FMS
subroutine mpp_error_mesg( routine, errormsg, errortype )
  character(len=*), intent(in) :: routine, errormsg
  integer,          intent(in) :: errortype

  call mpp_error( errortype, trim(routine)//': '//trim(errormsg) )
  return
end subroutine mpp_error_mesg

!#####################################################################
subroutine mpp_error_noargs()
  call mpp_error(FATAL)
end subroutine mpp_error_noargs

!#####################################################################
!--- makes a PE set out of a PE list. A PE list is an ordered list of PEs
!--- a PE set is a triad (start,log2stride,size) for SHMEM, an a communicator for MPI
!--- if stride is non-uniform or not a power of 2, 
!--- will return error (not required for MPI but enforced for uniformity)
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
 !synchronize PEs in list
subroutine mpp_sync( pelist )
  integer, intent(in), optional :: pelist(:)

  return
end subroutine mpp_sync

!#######################################################################
  !this is to check if current PE's outstanding puts are complete
  !but we can't use shmem_fence because we are actually waiting for
  !a remote PE to complete its get
subroutine mpp_sync_self( pelist )
  integer, intent(in), optional :: pelist(:)
  return
end subroutine mpp_sync_self

