#if defined(__sgi)
#define use_ftn_ieee_defs
#endif

module fpe_mod
!module to provide an uniform interface to individual (non-standard)
!libraries for controlling floating-point exception handling.

!this first version has a very simple interface with three routines:
!  logical function fpe_trap(flag)
!  subroutine fpe_disable(flag)
!  subroutine fpe_enable(flag)
!    integer, intent(in), optional :: flag
!in each case flag is an optional integer argument that is set to a sum of any
!of the public parameters:
!  FP_DIVIDE_BY_ZERO, FP_OPERAND_IS_NAN, FP_OVERFLOW, FP_UNDERFLOW
!an absent flag is equivalent to all of the above

!fpe_enable causes the program to abort if any exceptions in <flag> occur
!fpe_disable clears <flag>, so the program does not abort on those exceptions
!fpe_trap returns .TRUE. if any exceptions in <flag> have occurred
!  and the user retains control of the post-exception behaviour.
!  A side-effect of fpe_trap is to clear all exception bits.

!Interrupts will abort parallel executions (MPI or SMA) on a single system.

!This currently works on SGI systems only using the FTN_IEEE_DEFINITIONS
!compiler extension module. The module is also available on Cray systems
!but is not required there, as interrupts are enabled by default.

  use fms_mod, only: write_version_number

#ifdef use_ftn_ieee_defs
  use FTN_IEEE_DEFINITIONS
#endif
  implicit none
  private

  integer, parameter, public :: FP_DIVIDE_BY_ZERO=1, FP_OPERAND_IS_NAN=2, &
                                FP_OVERFLOW=4, FP_UNDERFLOW=8

  public :: fpe_init, fpe_trap, fpe_disable, fpe_enable

  logical :: module_is_initialized = .FALSE.

  character(len=128), private :: version= &
       '$Id: fpe.F90,v 2.3 2002/07/16 22:55:47 fms Exp $'
  character(len=128), private :: tagname= &
       '$Name: havana $'

  contains

    subroutine fpe_init ()

      if ( module_is_initialized ) return

      call write_version_number( version, tagname )

      module_is_initialized = .TRUE.

    end subroutine fpe_init

    function fpe_trap(flag)
!traps exceptions contained in the optional argument flag
      logical :: fpe_trap
      integer, intent(in), optional :: flag
      integer :: f
      logical :: status(0:3)

      if (.not.module_is_initialized ) call fpe_init ()

      f = FP_DIVIDE_BY_ZERO + FP_OPERAND_IS_NAN + FP_OVERFLOW + FP_UNDERFLOW
      if( PRESENT(flag) )f = flag
#ifdef use_ftn_ieee_defs
      status = TEST_IEEE_EXCEPTION( (/IEEE_XPTN_DIV_BY_ZERO, &
                                      IEEE_XPTN_INVALID_OPR, &
                                      IEEE_XPTN_OVERFLOW, &
                                      IEEE_XPTN_UNDERFLOW/) )
#else
      status = .FALSE.
#endif

      fpe_trap = .FALSE.
      if( BTEST(f,0) )fpe_trap = fpe_trap .OR. status(0)
      if( BTEST(f,1) )fpe_trap = fpe_trap .OR. status(1)
      if( BTEST(f,2) )fpe_trap = fpe_trap .OR. status(2)
      if( BTEST(f,3) )fpe_trap = fpe_trap .OR. status(3)

#ifdef use_ftn_ieee_defs
      call CLEAR_IEEE_EXCEPTION( (/IEEE_XPTN_DIV_BY_ZERO, &
                                   IEEE_XPTN_INVALID_OPR, &
                                   IEEE_XPTN_OVERFLOW, &
                                   IEEE_XPTN_UNDERFLOW/) )
#endif
      return
    end function fpe_trap

    subroutine fpe_disable(flag)
      integer, intent(in), optional :: flag
      integer :: f

      if (.not.module_is_initialized ) call fpe_init ()

      f = FP_DIVIDE_BY_ZERO + FP_OPERAND_IS_NAN + FP_OVERFLOW + FP_UNDERFLOW
      if( PRESENT(flag) )f = flag
#ifdef use_ftn_ieee_defs
      if( BTEST(f,0) )call DISABLE_IEEE_INTERRUPT(IEEE_NTPT_DIV_BY_ZERO)
      if( BTEST(f,1) )call DISABLE_IEEE_INTERRUPT(IEEE_NTPT_INVALID_OPR)
      if( BTEST(f,2) )call DISABLE_IEEE_INTERRUPT(IEEE_NTPT_OVERFLOW)
      if( BTEST(f,3) )call DISABLE_IEEE_INTERRUPT(IEEE_NTPT_UNDERFLOW)
#endif
      return
    end subroutine fpe_disable

    subroutine fpe_enable(flag)
      integer, intent(in), optional :: flag
      integer :: f

      if (.not.module_is_initialized ) call fpe_init ()

      f = FP_DIVIDE_BY_ZERO + FP_OPERAND_IS_NAN + FP_OVERFLOW + FP_UNDERFLOW
      if( PRESENT(flag) )f = flag
#ifdef use_ftn_ieee_defs
      if( BTEST(f,0) )call ENABLE_IEEE_INTERRUPT(IEEE_NTPT_DIV_BY_ZERO)
      if( BTEST(f,1) )call ENABLE_IEEE_INTERRUPT(IEEE_NTPT_INVALID_OPR)
      if( BTEST(f,2) )call ENABLE_IEEE_INTERRUPT(IEEE_NTPT_OVERFLOW)
      if( BTEST(f,3) )call ENABLE_IEEE_INTERRUPT(IEEE_NTPT_UNDERFLOW)
#endif
      return
    end subroutine fpe_enable
end module fpe_mod

#ifdef test_fpe
program fpetest
  use fpe_mod
  implicit none

  real :: a, b

  a = 1.
  b = 1.
  print *
  print *, 'Divide by zero with interrupts disabled...'
  call fpe_disable(FP_DIVIDE_BY_ZERO)
  print *, 'fpe_trap(FP_DIVIDE_BY_ZERO)=', fpe_trap(FP_DIVIDE_BY_ZERO)
  print *, 'divide by zero=', a/(a-b)
  print *, 'fpe_trap(FP_DIVIDE_BY_ZERO)=', fpe_trap(FP_DIVIDE_BY_ZERO)

  print *
  print *, 'Divide by zero with interrupts enabled...'
  print *, 'The next line should be an abort message:'
  call fpe_enable(FP_DIVIDE_BY_ZERO)
  print *, 'divide by zero=', a/(a-b)
  
end program fpetest
#endif

