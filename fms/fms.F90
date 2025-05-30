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
!> @defgroup fms_mod fms_mod
!> @ingroup fms
!! @brief The fms module provides routines that are commonly used
!!   by most FMS modules.
!> @author Bruce Wyman
!!
!> Here is a summary of the functions performed by routines
!!     in the fms module.
!!
!! 1. Output module version numbers to a common (<TT>log</TT>) file
!!     using a common format.<BR/>
!! 2. Open specific types of files common to many FMS modules.
!!     These include namelist files, restart files, and 32-bit IEEE
!!     data files. There also is a matching interface to close the files.
!!     If other file types are needed the <TT>mpp_open</TT> and <TT>mpp_close</TT>
!!     interfaces in module @ref mpp_io_mod must be used.<BR/>
!! 3. Read and write distributed data to simple native unformatted files.
!!     This type of file (called a restart file) is used to checkpoint
!!     model integrations for a subsequent restart of the run.<BR/>
!! 4. For convenience there are several routines published from
!!     the @ref mpp module. These are routines for getting processor
!!     numbers, commonly used I/O unit numbers, error handling, and timing sections of code.

!> @addtogroup fms_mod
!> @{
module fms_mod

!-----------------------------------------------------------------------
!
!         A collection of commonly used routines.
!
!  The routines are primarily I/O related, however, there also
!  exists several simple miscellaneous utility routines.
!
!-----------------------------------------------------------------------
!
!  file_exist         Checks the existence of the given file name.
!
!  check_nml_error    Checks the iostat argument that is returned after
!                     reading a namelist and determines if the error
!                     code is valid.
!
!  write_version_number  Prints to the log file (or a specified unit)
!                        the (cvs) version id string and (cvs) tag name.
!
!  error_mesg          Print notes, warnings and error messages,
!                      terminates program for error messages.
!                      (use error levels NOTE,WARNING,FATAL)
!
!  open_namelist_file  Opens namelist file for reading only.
!
!  open_restart_file   Opens a file that will be used for reading or writing
!                      restart files with native unformatted data.
!
!  open_ieee32_file    Opens a file that will be used for reading or writing
!                      unformatted 32-bit ieee data.
!
!  close_file          Closes a file that was opened using
!                      open_namelist_file, open_restart_file, or
!                      open_ieee32_file.
!
!  set_domain          Call this routine to internally store in fms_mod the
!                      domain2d data type prior to calling the distributed
!                      data I/O routines read_data and write_data.
!
!  read_data           Reads distributed data from a single threaded file.
!
!  write_data          Writes distributed data to a single threaded file.
!
!  fms_init            Initializes the fms module and also the
!                      mpp_io module (which initializes all mpp mods).
!                      Will be called automatically if the user does
!                      not call it.
!
!  fms_end             Calls mpp exit routines.
!
!  lowercase           Convert character strings to all lower case
!
!  uppercase           Convert character strings to all upper case
!
!  monotonic_array     Determines if the real input array has strictly
!                      monotonically increasing or decreasing values.
!
!  string_array_index  Match the input character string to a string
!                      in an array/list of character strings.
!
!-----------------------------------------------------------------------
!---- published routines from mpp_mod ----
!
!   mpp_error, NOTE, WARNING, FATAL
!   mpp_error_state
!   mpp_pe, mpp_npes, mpp_root_pe
!   stdin, stdout, stderr, stdlog
!   mpp_chksum
!
!   mpp_clock_id, mpp_clock_begin , mpp_clock_end
!   MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED
!   CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER,
!   CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA
!
!-----------------------------------------------------------------------

use          mpp_mod, only:  mpp_error, NOTE, WARNING, FATAL,    &
                             mpp_set_warn_level,                 &
                             mpp_transmit, ALL_PES,              &
                             mpp_pe, mpp_npes, mpp_root_pe,      &
                             mpp_sync, mpp_chksum,               &
                             mpp_clock_begin, mpp_clock_end,     &
                             mpp_clock_id, mpp_init, mpp_exit,   &
                             MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED, &
                             CLOCK_COMPONENT, CLOCK_SUBCOMPONENT,&
                             CLOCK_MODULE_DRIVER, CLOCK_MODULE,  &
                             CLOCK_ROUTINE, CLOCK_LOOP,          &
                             CLOCK_INFRA, mpp_clock_set_grain,   &
                             mpp_set_stack_size,                 &
                             stdin, stdout, stderr, stdlog,      &
                             mpp_error_state, lowercase,         &
                             uppercase, mpp_broadcast, input_nml_file, &
                             get_unit, read_input_nml

use  mpp_domains_mod, only:  domain2D, mpp_define_domains, &
                             mpp_update_domains, GLOBAL_DATA_DOMAIN, &
                             mpp_domains_init, mpp_domains_exit,     &
                             mpp_global_field, mpp_domains_set_stack_size,  &
                             mpp_get_compute_domain, mpp_get_global_domain, &
                             mpp_get_data_domain

use fms2_io_mod, only: fms2_io_init
use memutils_mod, only: print_memuse_stats, memutils_init
use grid2_mod, only: grid_init, grid_end
use fms_string_utils_mod, only: fms_c2f_string, fms_cstring2cpointer, string
use platform_mod, only: r4_kind, r8_kind

use, intrinsic :: iso_c_binding

implicit none
private

! routines for initialization and termination of module
public :: fms_init, fms_end

public ::check_nml_error, error_mesg, fms_error_handler

! version logging routine (originally from fms_io)
public :: write_version_number

! miscellaneous utilities (non i/o)
public :: lowercase, uppercase,        &
          string_array_index, monotonic_array

! public mpp interfaces
public :: mpp_error, NOTE, WARNING, FATAL, &
          mpp_error_state,                 &
          mpp_pe, mpp_npes, mpp_root_pe,   &
          stdin, stdout, stderr, stdlog,   &
          mpp_chksum, get_unit, read_input_nml
public :: input_nml_file
public :: mpp_clock_id, mpp_clock_begin, mpp_clock_end
public :: MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED
public :: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, &
          CLOCK_MODULE_DRIVER, CLOCK_MODULE,   &
          CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA
public :: fms_c2f_string, fms_cstring2cpointer
!public from the old fms_io but not exists here
public :: string

interface monotonic_array
  module procedure :: monotonic_array_r4, monotonic_array_r8
end interface monotonic_array

!Balaji
!this is published by fms and applied to any initialized clocks
!of course you can go and set the flag to SYNC or DETAILED by hand
integer, public :: clock_flag_default
!> @}
  !> Namelist read error values
  !> @ingroup fms_mod
  TYPE nml_errors_type
     INTEGER :: multipleNMLSinFile
     INTEGER :: badType1
     INTEGER :: badType2
     INTEGER :: missingVar
     INTEGER :: NotInFile
  END TYPE nml_errors_type
  TYPE(nml_errors_type), SAVE :: nml_errors
!> @addtogroup fms_mod
!> @{

!------ namelist interface -------
!------ adjustable severity level for warnings ------

  logical           :: read_all_pe   = .true. !< Read global data on all processors extracting local
                       !! part needed (TRUE) or read global data on PE0 and broadcast to all
                       !! PEs(FALSE).
  character(len=16) :: clock_grain = 'NONE' !< The level of clock granularity used for performance
                       !! timing sections of code. Possible values in order of increasing detail
                       !! are: 'NONE', 'COMPONENT', 'SUBCOMPONENT', 'MODULE_DRIVER', 'MODULE',
                       !! 'ROUTINE', 'LOOP', and 'INFRA'.  Code sections are defined using routines
                       !! in MPP module: mpp_clock_id, mpp_clock_begin, and mpp_clock_end. The fms
                       !! module makes these routines public. A list of timed code sections will be
                       !! printed to STDOUT. See the @ref mpp_mod module for more details.
  character(len=16) :: clock_flags='NONE' !< Possible values are 'NONE', 'SYNC', or 'DETAILED'.
                       !! SYNC will give accurate information on load balance of the clocked
                       !! portion of code. DETAILED also turns on detailed message-passing
                       !! performance diagnosis. Both SYNC and DETAILED will  work correctly on
                       !! innermost clock nest and distort outer clocks, and possibly the overall
                       !! code time. See the @ref mpp_mod module for more details.
  character(len=8)  :: warning_level = 'warning' !< Sets the termination condition for the WARNING
                       !! flag to interfaces error_mesg/mpp_error. set warning_level = 'fatal'
                       !! (program crashes for warning messages) or 'warning' (prints warning
                       !! message and continues).
  integer           :: stack_size = 0 !< The size in words of the MPP user stack. If stack_size > 0,
                       !! the following MPP routine is called: call mpp_set_stack_size (stack_size).
                       !! If stack_size = 0 (default) then the default size set by mpp_mod is used.
  integer           :: domains_stack_size = 0 !< The size in words of the MPP_DOMAINS user stack. If
                       !! domains_stack_size > 0, the following MPP_DOMAINS routine is called:
                       !! call mpp_domains_set_stack_size (domains_stack_size). If
                       !! domains_stack_size = 0 (default) then the default size set by
                       !! @ref mpp_domains_mod is used.
  logical, public   :: print_memory_usage = .FALSE. !< If set to .TRUE., memory usage statistics
                       !! will be printed at various points in the code. It is used to study memory
                       !! usage, e.g to detect memory leaks.

!------ namelist interface -------

  namelist /fms_nml/  read_all_pe, clock_grain, clock_flags,         &
                      warning_level, stack_size, domains_stack_size, &
                      print_memory_usage

!   ---- private data for check_nml_error ----

   integer, private :: num_nml_error_codes, nml_error_codes(20)
   logical, private :: do_nml_error_init = .true.
   private  nml_error_init


!  ---- version number -----

! Include variable "version" to be written to log file.
#include<file_version.h>

  logical :: module_is_initialized = .FALSE.

!> @}

!> @addtogroup fms_mod
!> @{
contains

!#######################################################################

!> @brief Initializes the FMS module and also calls the initialization routines for all
!!     modules in the MPP package. Will be called automatically if the user does
!!     not call it.
!! @details Initialization routine for the fms module. It also calls initialization routines
!!      for the mpp, mpp_domains, and mpp_io modules. Although this routine
!!      will be called automatically by other fms_mod routines, users should
!!      explicitly call fms_init. If this routine is called more than once it will
!!      return silently. There are no arguments.
!!
!> @throws FATAL, invalid entry for namelist variable warning_level
!! The namelist variable warning_level must be either 'fatal' or 'warning'(case-insensitive)
!!
!> @throws FATAL, invalid entry for namelist variable clock_grain
!! The namelist variable clock_grain must be one of the following values:
!! 'NONE', 'COMPONENT', 'SUBCOMPONENT', 'MODULE_DRIVER', 'MODULE', 'ROUTINE',
!! 'LOOP', or 'INFRA' (case-insensitive).
subroutine fms_init (localcomm, alt_input_nml_path)

!--- needed to output the version number of constants_mod to the logfile ---
 use constants_mod, only: constants_version=>version !pjp: PI not computed
 interface
    subroutine maximize_system_stacksize_limit() bind(C)
    end subroutine
 end interface

 integer, intent(in), optional :: localcomm
 character(len=*), intent(in), optional :: alt_input_nml_path
 integer :: ierr, io
 integer :: logunitnum
 integer :: stdout_unit !< Unit number for the stdout file

    if (module_is_initialized) return    ! return silently if already called
    module_is_initialized = .true.

!---- Raise the system stack size limit to its maximum permissible value ----
    call maximize_system_stacksize_limit

!---- initialize mpp routines ----
    if(present(localcomm)) then
       if(present(alt_input_nml_path)) then
          call mpp_init(localcomm=localcomm, alt_input_nml_path=alt_input_nml_path)
       else
          call mpp_init(localcomm=localcomm)
       endif
    else
       if(present(alt_input_nml_path)) then
          call mpp_init(alt_input_nml_path=alt_input_nml_path)
       else
          call mpp_init()
       endif
    endif
    call mpp_domains_init()
    call fms2_io_init()
#ifdef use_deprecated_io
      call mpp_error(NOTE, "fms_io_init: fms_io HAS BEEN DEPRECATED! "//&
                           "PLEASE REMOVE -Duse_deprecated_io FROM YOUR COMPILE FLAGS "// &
                           "AND MOVE TO FMS2_IO. CONTACT YOUR MODEL LIASISON IF YOU NEED "// &
                           "ASSISTANCE")
#endif
    logunitnum = stdlog()
!---- read namelist input ----

    call nml_error_init()  ! first initialize namelist iostat error codes

    read (input_nml_file, fms_nml, iostat=io)
    ierr = check_nml_error(io,'fms_nml')

!---- define mpp stack sizes if non-zero -----

    if (        stack_size > 0) call         mpp_set_stack_size (        stack_size)
    if (domains_stack_size > 0) call mpp_domains_set_stack_size (domains_stack_size)

!---- set severity level for warnings ----

    select case( trim(lowercase(warning_level)) )
    case( 'fatal' )
        call mpp_set_warn_level ( FATAL )
    case( 'warning' )
        call mpp_set_warn_level ( WARNING )
    case default
        call error_mesg ( 'fms_init',  &
             'invalid entry for namelist variable warning_level', FATAL )
    end select

!--- set granularity for timing code sections ---

    select case( trim(uppercase(clock_grain)) )
    case( 'NONE' )
        call mpp_clock_set_grain (0)
    case( 'COMPONENT' )
        call mpp_clock_set_grain (CLOCK_COMPONENT)
    case( 'SUBCOMPONENT' )
        call mpp_clock_set_grain (CLOCK_SUBCOMPONENT)
    case( 'MODULE_DRIVER' )
        call mpp_clock_set_grain (CLOCK_MODULE_DRIVER)
    case( 'MODULE' )
        call mpp_clock_set_grain (CLOCK_MODULE)
    case( 'ROUTINE' )
        call mpp_clock_set_grain (CLOCK_ROUTINE)
    case( 'LOOP' )
        call mpp_clock_set_grain (CLOCK_LOOP)
    case( 'INFRA' )
        call mpp_clock_set_grain (CLOCK_INFRA)
    case default
        call error_mesg ( 'fms_init',  &
             'invalid entry for namelist variable clock_grain', FATAL )
    end select
!Balaji
    select case( trim(uppercase(clock_flags)) )
    case( 'NONE' )
       clock_flag_default = 0
    case( 'SYNC' )
       clock_flag_default = MPP_CLOCK_SYNC
    case( 'DETAILED' )
       clock_flag_default = MPP_CLOCK_DETAILED
    case default
       call error_mesg ( 'fms_init',  &
            'invalid entry for namelist variable clock_flags', FATAL )
   end select

!--- write version info and namelist to logfile ---

    call write_version_number("FMS_MOD", version)
    if (mpp_pe() == mpp_root_pe()) then
      stdout_unit = stdlog()
      write (stdout_unit, nml=fms_nml)
      write (stdout_unit,*) 'nml_error_codes=', nml_error_codes(1:num_nml_error_codes)
    endif

    call memutils_init( print_memory_usage )
    call print_memuse_stats('fms_init')

!--- output version information constants to the logfile
    call write_version_number("CONSTANTS_MOD", constants_version)
    call grid_init

end subroutine fms_init

!#######################################################################

!> @brief Calls the termination routines for all modules in the MPP package.
!!
!> Termination routine for the fms module. It also calls destructor routines
!! for the mpp, mpp_domains, and mpp_io modules. If this routine is called
!! more than once it will return silently. There are no arguments.
subroutine fms_end ( )

    if (.not.module_is_initialized) return  ! return silently
!    call fms_io_exit  ! now called from coupler_end
    call grid_end
    call mpp_domains_exit
    call mpp_exit
    module_is_initialized =.FALSE.

end subroutine fms_end

!#######################################################################

 !> @brief Print notes, warnings and error messages; terminates program for warning
 !! and error messages. Usage of @ref mpp_error is preferable. (use error levels NOTE,WARNING,FATAL, see example below)
 !! @details Print notes, warnings and error messages; and terminates the program for
 !!     error messages. This routine is a wrapper around mpp_error, and is provided
 !!     for backward compatibility. This module also publishes mpp_error,
 !!      <B>users should try to use the mpp_error interface</B>.
 !!
 !! <br>Example usage:
 !! @code{.F90}
 !! use fms_mod, only: error_mesg, FATAL, NOTE
 !! call error_mesg ('fms_mod', 'initialization not called', FATAL)
 !! call error_mesg ('fms_mod', 'fms_mod message', NOTE)
 !! @endcode
 subroutine error_mesg (routine, message, level)
  character(len=*), intent(in) :: routine !< Routine name where the warning or error has occurred.
  character(len=*), intent(in) :: message !< Warning or error message to be printed.
  integer,          intent(in) :: level !< Level of severity; set to NOTE, WARNING, or FATAL Termination always occurs
                                        !! for FATAL, never for NOTE, and is settable for WARNING (see namelist).

!  input:
!      routine   name of the calling routine (character string)
!      message   message written to output   (character string)
!      level     set to NOTE, MESSAGE, or FATAL (integer)

    if (.not.module_is_initialized) call fms_init ( )
    call mpp_error ( routine, message, level )

 end subroutine error_mesg

!#######################################################################

 !> @brief Facilitates the control of fatal error conditions
 !! @details When err_msg is present, message is copied into err_msg
 !!     and the function returns a value of .true.
 !!     Otherwise calls mpp_error to terminate execution.
 !!     The intended use is as shown below.
 !! @returns true when err_msg is present
 !! @code{.F90}
 !! if(fms_error_handler(routine, message, err_msg)) return
 !! @endcode
 function fms_error_handler(routine, message, err_msg)

 logical :: fms_error_handler
 character(len=*), intent(in) :: routine !< Routine name where the fatal error has occurred.
 character(len=*), intent(in) :: message !< fatal error message to be printed.
 character(len=*), intent(out), optional :: err_msg !< When err_msg is present: err_msg = message

 fms_error_handler = .false.
 if(present(err_msg)) then
   err_msg = message
   fms_error_handler = .true.
 else
   call mpp_error(trim(routine),trim(message),FATAL)
 endif

 end function fms_error_handler

! used to check the iostat argument that is
! returned after reading a namelist
! see the online documentation for how this routine might be used

  !> @brief Checks the iostat argument that is returned after reading a namelist
  !!     and determines if the error code is valid.
  !! @return This function returns the input iostat value (integer) if it is an
  !!     allowable error code. If the iostat error code is not
  !!     allowable, an error message is printed and the program terminated.
  !! @details The FMS allows multiple namelist records to reside in the same file.
  !!     Use this interface to check the iostat argument that is returned after
  !!     reading a record from the namelist file. If an invalid iostat value
  !!     is detected this routine will produce a fatal error. See the NOTE below.
  !!
  !!     Some compilers will return non-zero iostat values when reading through
  !!     files with multiple namelist. This routine
  !!     will try skip these errors and only terminate for true namelist errors.
  !!
  !!     <br>Examples<br>
  !!
  !!       The following example checks if a file exists, reads a namelist input
  !!       from that file, and checks for errors in that
  !!       namelist. When the correct namelist is read and it has no errors the
  !!       routine check_nml_error will return zero and the while loop will exit.
  !!       This code segment should be used to read namelist files.
  !!       @code{.F90}
  !!         integer :: ierr, io
  !!
  !!         read (input_nml_file, fms_nml, iostat=io)
  !!         ierr = check_nml_error(io,'fms_nml')
  !!       @endcode
  !! @throws FATAL, Unknown error while reading namelist ...., (IOSTAT = ####)
  !! There was an error reading the namelist specified. Carefully examine all namelist and variables
  !! for anything incorrect (e.g. malformed, hidden characters).
  !!
  !! @throws FATAL, Unknown namelist, or mistyped namelist variable in namelist ...., (IOSTAT = ####)
  !! The name list given doesn't exist in the namelist file, or a variable in the namelist is
  !! mistyped or isn't a namelist variable.
  INTEGER FUNCTION check_nml_error(IOSTAT, NML_NAME)
    INTEGER, INTENT(in) :: IOSTAT !< The iostat value returned when reading a namelist record.
    CHARACTER(len=*), INTENT(in) :: NML_NAME !< The name of the namelist. This name will be printed if an error is
                                             !! encountered, otherwise the name is not used.

    CHARACTER(len=256) :: err_str

    IF ( .NOT.module_is_initialized) CALL fms_init()

    check_nml_error = IOSTAT

    ! Return on valid IOSTAT values
    IF ( IOSTAT <= 0 .OR.&
       & IOSTAT == nml_errors%multipleNMLSinFile .OR.&
       & IOSTAT == nml_errors%NotInFile) RETURN

    ! Everything else is a FATAL
    IF ( (IOSTAT == nml_errors%badType1 .OR. IOSTAT == nml_errors%badType2) .OR. IOSTAT == nml_errors%missingVar ) THEN
       WRITE (err_str,*) 'Unknown namelist, or mistyped namelist variable in namelist ',TRIM(NML_NAME),', &
             &  (IOSTAT = ',IOSTAT,')'
       CALL error_mesg ('check_nml_error in fms_mod', err_str, FATAL)
       CALL mpp_sync()
    ELSE
       WRITE (err_str,*) 'Unknown error while reading namelist ',TRIM(NML_NAME),', (IOSTAT = ',IOSTAT,')'
       CALL error_mesg ('check_nml_error in fms_mod', err_str, FATAL)
       CALL mpp_sync()
    END IF
  END FUNCTION check_nml_error

!-----------------------------------------------------------------------
!   private routine for initializing allowable error codes

  !> @brief Determines the IOSTAT error value for some common Namelist errors.
  !!   Also checks if the compiler returns a non-zero status if there are
  !!   multiple namelist records in a single file.
  SUBROUTINE nml_error_init
    ! Determines the IOSTAT error value for some common Namelist errors.
    ! Also checks if the compiler returns a non-zero status if there are
    ! multiple namelist records in a single file.
    INTEGER, PARAMETER :: unit_begin = 20, unit_end = 1024
    INTEGER :: fileunit, io_stat
    INTEGER, DIMENSION(5) :: nml_iostats
    LOGICAL :: opened

    ! Variables for sample namelists
    INTEGER :: i1 !< Variables for sample namelists
    INTEGER :: i2 !< Variables for sample namelists
    REAL :: r1, r2
    LOGICAL :: l1
    NAMELIST /a_nml/ i1, r1
    NAMELIST /b_nml/ i2, r2, l1
    NAMELIST /badType1_nml/ i1, r1
    NAMELIST /badType2_nml/ i1, r1
    NAMELIST /missingVar_nml/ i2, r2
    NAMELIST /not_in_file_nml/ i2, r2

    ! Initialize the sample namelist variables
    i1 = 1
    i2 = 2
    r1 = 1.0
    r2 = 2.0
    l1 = .FALSE.

    ! Create a dummy namelist file
    IF ( mpp_pe() == mpp_root_pe() ) THEN
       ! Find a free file unit for a scratch file
       file_opened: DO fileunit = unit_begin, unit_end
          INQUIRE(UNIT=fileunit, OPENED=opened)
          IF ( .NOT.opened ) EXIT file_opened
       END DO file_opened

#if defined(__PGI) || defined(_CRAYFTN)
       OPEN (UNIT=fileunit, FILE='_read_error.nml', IOSTAT=io_stat)
#else
       OPEN (UNIT=fileunit, STATUS='SCRATCH', IOSTAT=io_stat)
#endif

       ! Write sample namelist to the SCRATCH file.
       WRITE (UNIT=fileunit, NML=a_nml, IOSTAT=io_stat)
       WRITE (UNIT=fileunit, NML=b_nml, IOSTAT=io_stat)
       WRITE (UNIT=fileunit, IOSTAT=io_stat, FMT='(/,"&badType1_nml  i1=1, r1=''bad'' /",/)')
       WRITE (UNIT=fileunit, IOSTAT=io_stat, FMT='(/,"&badType2_nml  i1=1, r1=.true. /",/)')
       WRITE (UNIT=fileunit, IOSTAT=io_stat, FMT='(/,"&missingVar_nml  i2=1, r2=1.0e0, l1=.true. /",/)')

       ! Rewind for reading
       REWIND(UNIT=fileunit)

       ! Read the second namelist from the file -- check for namelist bug
       READ (UNIT=fileunit, NML=b_nml, IOSTAT=nml_iostats(1))
       REWIND(UNIT=fileunit)

       ! Read in bad type 1 --- Some compilers treat the string cast differently
       READ (UNIT=fileunit, NML=badType1_nml, IOSTAT=nml_iostats(2))
       REWIND(UNIT=fileunit)

       ! Read in bad type 2
       READ (UNIT=fileunit, NML=badType2_nml, IOSTAT=nml_iostats(3))
       REWIND(UNIT=fileunit)

       ! Read in missing variable/misstyped
       READ (UNIT=fileunit, NML=missingVar_nml, IOSTAT=nml_iostats(4))
       REWIND(UNIT=fileunit)

       ! Code for namelist not in file
       READ (UNIT=fileunit, NML=not_in_file_nml, IOSTAT=nml_iostats(5))

       ! Done, close file
       CLOSE (UNIT=fileunit)

       ! Some compilers don't handle the type casting as well as we would like.
       IF ( nml_iostats(2) * nml_iostats(3) .EQ. 0 ) THEN
          IF ( nml_iostats(2) .NE. 0 .AND. nml_iostats(3) .EQ. 0 ) THEN
             nml_iostats(3) = nml_iostats(2)
          ELSE IF ( nml_iostats(2) .EQ. 0 .AND. nml_iostats(3) .NE.0 ) THEN
             nml_iostats(2) = nml_iostats(3)
          ELSE
             nml_iostats(2) = nml_iostats(4)
             nml_iostats(2) = nml_iostats(4)
          END IF
       END IF
    END IF

    ! Broadcast nml_errors
    CALL mpp_broadcast(nml_iostats,5,mpp_root_pe())
    nml_errors%multipleNMLSinFile = nml_iostats(1)
    nml_errors%badType1 = nml_iostats(2)
    nml_errors%badType2 = nml_iostats(3)
    nml_errors%missingVar = nml_iostats(4)
    nml_errors%NotInFile = nml_iostats(5)

    do_nml_error_init = .FALSE.
  END SUBROUTINE nml_error_init

!#######################################################################

!> @brief match the input character string to a string
!!     in an array/list of character strings
!! @return If an exact match was found then true is returned, otherwise false is returned.
!! @details Tries to find a match for a character string in a list of character strings.
!!      The match is case sensitive and disregards blank characters to the right of
!!      the string.
!!
!!      <br>Examples<br>
!!      @code{.F90}
!!       string = "def"
!!       string_array = (/ "abcd", "def ", "fghi" /)
!!
!!       string_array_index ( string, string_array, index )
!!      @endcode
!!       Returns: TRUE, index = 2
function string_array_index ( string, string_array, index ) result (found)
character(len=*),  intent(in)  :: string !< Character string of arbitrary length.
character(len=*),  intent(in)  :: string_array(:) !< Array/list of character strings.
integer, optional, intent(out) :: index !< The index of string_array where the first match was found. If
                                        !! no match was found then index = 0.
logical :: found !< If an exact match was found then TRUE is returned, otherwise FALSE is returned.
integer :: i

! initialize this function to false
! loop thru string_array and exit when a match is found

  found = .false.
  if (present(index)) index = 0

  do i = 1, size(string_array(:))
    ! found a string match ?
    if ( trim(string) == trim(string_array(i)) ) then
         found = .true.
         if (present(index)) index = i
         exit
    endif
  enddo

end function string_array_index

!#######################################################################
!> @brief Prints to the log file (or a specified unit) the version id string and
!!  tag name.
subroutine write_version_number (version, tag, unit)
  character(len=*), intent(in) :: version !> string that contains routine name
  character(len=*), intent(in), optional :: tag !> tag name that code was checked out with
  integer,          intent(in), optional :: unit !> alternate unit number to direct output,
                                                 !! defaults to stdlog
  integer :: logunit

  if (.not.module_is_initialized) call fms_init ( )

  logunit = stdlog()

  if (present(unit)) then
    logunit = unit
  else
    ! only allow stdlog messages on root pe
    if ( mpp_pe() /= mpp_root_pe() ) return
  endif

  if (present(tag)) then
    write (logunit,'(/,80("="),/(a))') trim(version), trim(tag)
  else
    write (logunit,'(/,80("="),/(a))') trim(version)
  endif

end subroutine write_version_number

#include "fms_r4.fh"
#include "fms_r8.fh"

end module fms_mod
! <INFO>
!   <BUG>
!     Namelist error checking may not work correctly with some compilers.
!
!     Users should beware when mixing Fortran reads and read_data calls. If a
!     Fortran read follows read_data and namelist variable read_all_pe = FALSE
!     (not the default), then the code will fail. It is safest if Fortran reads
!     precede calls to read_data.
!   </BUG>
!   <ERROR MSG="unexpected EOF" STATUS="FATAL">
!     An unexpected end-of-file was encountered in a read_data call.
!     You may want to use the optional end argument to detect the EOF.
!   </ERROR>
!   <NOTE>
!     1) If the <B>MPP</B> or <B>MPP_DOMAINS</B> stack size is exceeded the
!     program will terminate after printing the required size.
!
!     2) When running on a very small number of processors or for high
!     resolution models the default domains_stack_size will
!     probably be insufficient.
!
!     3) The following performance routines in the <B>MPP</B> module are published by this module.
!<PRE>
!        mpp_clock_id, mpp_clock_begin, mpp_clock_end
!</PRE>
!        and associated parameters that are published:
!<PRE>
!        MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED, CLOCK_COMPONENT, CLOCK_SUBCOMPONENT,
!        CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA
!</PRE>
!
!     4) Here is an example of how to time a section of code.<BR/>
!<PRE>
!          use fms_mod, only: mpp_clock_id, mpp_clock_begin, &
!                             mpp_clock_end. MPP_CLOCK_SYNC, &
!                             CLOCK_MODULE_DRIVER
!          integer :: id_mycode
!
!          id_mycode = mpp_clock_id ('mycode loop', flags=MPP_CLOCK_SYNC, grain=CLOCK_MODULE_DRIVER)
!          call mpp_clock_begin (id_mycode)
!                        :
!                        :
!           ~~ this code will be timed ~~
!                        :
!                        :
!          call mpp_clock_end (id_mycode)
! </PRE>
!        Note: <TT>CLOCK_MODULE_DRIVER</TT> can be replaced with
!        CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE,
!        CLOCK_LOOP, or CLOCK_INFRA.
!
!   </NOTE>
!   <FUTURE>
!     NetCDF facilities for reading and writing restart files and (IEEE32)
!       data files.
!    </FUTURE>
!    <FUTURE>
!     May possible split the FMS module into two modules.
!
!      i.general utilities (FMS_MOD) <BR/>
!     ii.I/O utilities (FMS_IO_MOD)
!    </FUTURE>
! </INFO>
!> @}
! close documentation grouping
