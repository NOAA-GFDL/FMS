#if defined(__sgi)
#define use_ftn_ieee_defs
#endif

module fpe_mod

! <CONTACT EMAIL="vb@gfdl.noaa.gov">
!    V. Balaji
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   The module <TT>fpe</TT> provides a set of simple calls for controlling
!   code behaviour under floating-point exceptions.
! </OVERVIEW>

! <DESCRIPTION>
!   Numbers on digital computers are represented as a series of bits,
!   and thus intrinsically only deal with integers. A <I>number model</I>
!   for "real" numbers is generally used to provide a representation of
!   non-integers. Since only a finite number of bits is used, this
!   representation does not actually map onto the real number set, but a
!   finite set of those numbers admitting a terminating representation
!   (i.e excluding irrational numbers and non-terminating fractions). This
!   number model is now almost universally based on <I>floating-point (FP)
!   numbers</I>, using a mantissa and an exponent. The IEEE has
!   specified a standard FP number model, which is used on most, but not
!   all, modern computers.
!   When FP computations produce numbers that cannot be represented in the
!   number model, an <I>exception</I> is said to have occurred. Exceptions
!   fall into several classes, of which the most common include:
!
!  <B>Inexact result</B>: FP result has no exact representation
!  and must be rounded. This occurs extremely frequently: in fact the
!  probability that an FP calculation has an exact result is in theory
!  vanishingly small.<BR/>
!  <B>Overflow</B>: Absolute value of FP result is larger than
!  the largest FP number representable. The F90 intrinsic <TT>HUGE()</TT>
!  returns this value.<BR/>
!  <B>Underflow</B>: Absolute value of FP result is smaller than
!  the smallest FP number representable. The F90 intrinsic
!  <TT>TINY()</TT> returns this value.<BR/>
!  <B>Divide by zero</B>: FP division has a 0 in the denominator.<BR/>
!  <B>Invalid operand</B>: FP operation is attempting to treat a
!  series of bits that is not a valid FP number. Invalid FP numbers are
!  generally called <I>NaN</I>s (not-a-number). This exception class may
!  also cover other invalid operands, such as a negative argument to a FP
!  square-root operation.
!  
!  Computing units maintain a word in memory in order to track the
!  occurrence of exceptions. Each bit in this word signals the occurrence
!  of a different class of exceptions, such as those listed above. This
!  word must be tested to know if an exception has occurred. The result
!  of this test can be used to control subsequent program behaviour.
!  Principally, we might wish to interrupt the program if certain classes
!  of exceptions occur, though other more refined consequences may also
!  be imagined.
!  
!  There is a computational burden associated with testing for
!  exceptions, and also one with continuing the execution of a program
!  after the computation has gone away. This module provides a standard
!  interface to a simple range of behaviour under exceptions. This
!  includes enabling and disabling interrupts for certain classes of
!  exceptions, and a trap function that merely signals to the user that
!  an exception has occurred.

! </DESCRIPTION>

! <PUBLIC>
!     The public interfaces to the <TT>fpe</TT> module are described here
!     in alphabetical order. The <TT>fpe</TT> module defines four classes of
!     exceptions: <TT>FP_DIVIDE_BY_ZERO, FP_OPERAND_IS_NAN, FP_OVERFLOW,
!     FP_UNDERFLOW</TT>, which are public integer parameters use-associated
!     from the module. The argument <TT>flag</TT> to any of the routines
!     below can be a sum of any of these. An absent <TT>flag</TT> argument
!     is equivalent to all of the exceptions, i.e equivalent to
!     <TT>flag=FP_DIVIDE_BY_ZERO+FP_OPERAND_IS_NAN+FP_OVERFLOW+FP_UNDERFLOW</TT>.
! </PUBLIC>

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
       '$Id: fpe.F90,v 2.4 2003/04/09 21:17:04 fms Exp $'
  character(len=128), private :: tagname= &
       '$Name: inchon $'

  contains

    subroutine fpe_init ()

      if ( module_is_initialized ) return

      call write_version_number( version, tagname )

      module_is_initialized = .TRUE.

    end subroutine fpe_init

! <FUNCTION NAME="fpe_trap">

!   <OVERVIEW>
!     Set trap for classes of exceptions. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     Traps exceptions contained in the optional argument flag.
!   </DESCRIPTION>
!   <TEMPLATE>
!     logical function fpe_trap(flag)
!   </TEMPLATE>

!   <IN NAME="flag" TYPE="integer"> </IN>
!   <OUT>
!     Returns <TT>.TRUE.</TT> if any of the exceptions
!     in <TT>flag</TT> has occurred. This may be an alternative to the
!     computational burden of <TT>fpe_enable</TT>, since the user retains
!     control of when to test for exceptions.
!   </OUT>

    function fpe_trap(flag)
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
! </FUNCTION>

! <SUBROUTINE NAME="fpe_disable">

!   <OVERVIEW>
!     Disable interrupts for classes of exceptions.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Causes the program to continue execution even if one of the exceptions
!     in flag has occurred, until a subsequent call to fpe_enable.
!   </DESCRIPTION>
!   <TEMPLATE> call fpe_disable(flag)</TEMPLATE>
!   <IN NAME="flag" TYPE="integer"> </IN>

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
! </SUBROUTINE>

! <SUBROUTINE NAME="fpe_enable">
      
!   <OVERVIEW>
!     Enable interrupts for classes of exceptions.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Causes the program to abort if any of the exceptions in flag occurs.
!     This involves frequent exception testing and may slow down the program.
!   </DESCRIPTION>
!   <TEMPLATE> call fpe_enable(flag)</TEMPLATE>
!   <IN NAME="flag" TYPE="integer"> </IN>

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
! </SUBROUTINE>

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

! <INFO>

!   <COMPILER NAME="">     
!     <h4>COMPILING AND LINKING SOURCE</h4>
!     Any module or program unit using the <TT>fpe</TT> module must
!     contain the line 
!     
!     <PRE>
!     use fpe
!     </PRE>
!     
!     The source file for the <TT>fpe</TT> module is <LINK SRC
!     ="ftp://ftp.gfdl.noaa.gov/pub/vb/utils/fpe.F90">fpe.F90</LINK>.
!     
!     Compiling with the cpp flag <TT>test_fpe</TT> turned on:
!     
!     <PRE>
!     f90 -Dtest_fpe fpe.F90
!     </PRE>
!     
!     will produce a program that will exercise certain portions of the
!     the <TT>fpe</TT>  module.

!   </COMPILER>
!   <PRECOMP FLAG="">      
!     The <TT>fpe</TT> module was written specifically for the SGIs using
!     the MIPSpro f90 compiler. They use the non-standard
!     <TT>FTN_IEEE_DEFINITIONS</TT> module that comes with this compiler and
!     provides exception handling functionality following the IEEE FP number
!     model.
!     
!     This will shortly be extended to various compilers on the Linux
!     platform. The module is deemed unnecessary on Crays, on which GFDL
!     seems to have got by quite satisfactorily for many years without this
!     degree of user control.
!     
!     Tim Yeager (ty@gfdl.noaa.gov) has written an
!     excellent introduction to FPE behaviour for GFDL users, including
!     various run-time and compile-time options for controlling this
!     behaviour. Slides from his recent talk on the subject, and associated
!     files, are available in <TT>/home/ty/doc/trap.talk</TT>.

!   </PRECOMP> 
!   <LOADER FLAG="">       
!     GFDL users can copy the file <TT>/net/vb/public/utils/fpe.F90</TT>. 
!     External users can download the source <LINK SRC =
!     "ftp://ftp.gfdl.noaa.gov/pub/vb/utils/fpe.F90">here</LINK>. The current 
!     public version number is 2.2. 
!   </LOADER>
!   <TESTPROGRAM NAME="">  text </TESTPROGRAM>
!   <NOTE>
!     The <TT>fpe</TT> module defines four classes of exceptions:
!  <TT>FP_DIVIDE_BY_ZERO, FP_OPERAND_IS_NAN, FP_OVERFLOW,
!  FP_UNDERFLOW</TT>. These are <TT>integer</TT>s and may be summed to
!  designate multiple classes.
!  
!  The calls <TT>fpe_enable</TT> and <TT>fpe_disable</TT> may be used
!  to turn on FP error interrupts for different code
!  sections. <TT>fpe_trap</TT> returns <TT>.TRUE.</TT> if an exception
!  has occurred. Consider the test program:
!  
!  <PRE>
!  program fpetest
!    use fpe
!    implicit none
!  
!    real :: a, b
!  
!    a = 1.
!    b = 1.
!    print *
!    print *, 'Divide by zero with interrupts disabled...'
!    call fpe_disable(FP_DIVIDE_BY_ZERO)
!    print *, 'fpe_trap(FP_DIVIDE_BY_ZERO)=', fpe_trap(FP_DIVIDE_BY_ZERO)
!    print *, 'divide by zero=', a/(a-b)
!    print *, 'fpe_trap(FP_DIVIDE_BY_ZERO)=', fpe_trap(FP_DIVIDE_BY_ZERO)
!  
!    print *
!    print *, 'Divide by zero with interrupts enabled...'
!    print *, 'The next line should be an abort message:'
!    call fpe_enable(FP_DIVIDE_BY_ZERO)
!    print *, 'divide by zero=', a/(a-b)
!  
!  end program fpetest
!  </PRE>
!  
!  This produces (on SGIs) the output:
!  
!  <PRE>
!   Divide by zero with interrupts disabled...
!   fpe_trap(FP_DIVIDE_BY_ZERO)= F
!   divide by zero= Infinity
!   fpe_trap(FP_DIVIDE_BY_ZERO)= T
!  
!   Divide by zero with interrupts enabled...
!   The next line should be an abort message:
!   Floating Exception
!   Abort
!  </PRE>
!  
!  This example provides a comprehensive test case for the FPE module.
!  
!  In the first computation of <TT>1./(1.-1.)</TT>, FPE interrupts
!  have been disabled.
!  Before the computation the <TT>FP_DIVIDE_BY_ZERO</TT> exception is
!  shown as not having occurred. Subsequently this exception bit has been
!  set to <TT>.TRUE.</TT> In the second instance, FPE interrupts have
!  been enabled, and cause the program to abort.
!   </NOTE>

! </INFO>
