
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
!  monotonic_array     Determines if the real input array has
!                      monotonically increasing or decreasing values.
!
!  string_array_index  Match the input character string to a string
!                      in an array/list of character strings.
!
!  mpp_clock_init      Sets up a identifier for performance timing
!                        (similar to mpp_clock_id)
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
!   mpp_clock_begin , mpp_clock_end
!   MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED
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
                             mpp_set_stack_size,                 &
                             stdin, stdout, stderr, stdlog,      &
                             mpp_error_state

use  mpp_domains_mod, only:  domain2D, mpp_define_domains, &
                             mpp_update_domains, GLOBAL_DATA_DOMAIN, &
                             mpp_domains_init, mpp_domains_exit,     &
                             mpp_global_field, mpp_domains_set_stack_size,  &
                             mpp_get_compute_domain, mpp_get_global_domain, &
                             mpp_get_data_domain

use       mpp_io_mod, only:  mpp_io_init, mpp_open, mpp_close,         &
                       MPP_ASCII, MPP_NATIVE, MPP_IEEE32, MPP_NETCDF,  &
                       MPP_RDONLY, MPP_WRONLY, MPP_APPEND, MPP_OVERWR, &
                       MPP_SEQUENTIAL, MPP_DIRECT,                     &
                       MPP_SINGLE, MPP_MULTI, MPP_DELETE, mpp_io_exit

implicit none
private

! routines for initialization and termination of module
public :: fms_init, fms_end

! routines for opening/closing specific types of file
public :: open_namelist_file, open_restart_file, &
          open_ieee32_file, close_file

! routines for reading/writing distributed data
public :: set_domain, read_data, write_data
public :: get_domain_decomp

! miscellaneous i/o routines
public :: file_exist, check_nml_error,      &
          write_version_number, error_mesg

! miscellaneous utilities (non i/o)
public :: lowercase, uppercase,                &
          string_array_index, monotonic_array, &
          mpp_clock_init

! public mpp interfaces
public :: mpp_error, NOTE, WARNING, FATAL, &
          mpp_error_state,                 &
          mpp_pe, mpp_npes, mpp_root_pe,   &
          stdin, stdout, stderr, stdlog,   &
          mpp_chksum
public :: mpp_clock_begin, mpp_clock_end
public :: MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED
           

! overloaded distributed data i/o routines

interface read_data
  module procedure read_data_2d, read_ldata_2d, read_idata_2d
  module procedure read_data_3d, read_data_4d
  module procedure read_cdata_2d,read_cdata_3d,read_cdata_4d
end interface

interface write_data
  module procedure write_data_2d, write_ldata_2d, write_idata_2d
  module procedure write_data_3d, write_data_4d
  module procedure write_cdata_2d,write_cdata_3d,write_cdata_4d
end interface

!------ namelist interface -------
!------ adjustable severity level for warnings ------

  integer           :: timing_level  = 0
  logical           :: read_all_pe   = .true.
  character(len=8)  :: warning_level = 'warning'
  character(len=64) :: iospec_ieee32 = '-N ieee_32'
  integer           :: stack_size = 0
  integer           :: domains_stack_size = 0

  namelist /fms_nml/  timing_level, read_all_pe,    &
                      warning_level, iospec_ieee32, &
                      stack_size, domains_stack_size


!   ---- private data for check_nml_error ----

   integer, private :: num_nml_error_codes, nml_error_codes(20)
   logical, private :: do_nml_error_init = .true.
   private  nml_error_init


!------ private data, pointer to current 2d domain ------

  type(domain2D), pointer :: Domain

  integer :: is,ie,js,je      ! compute domain
  integer :: isd,ied,jsd,jed  ! data domain
  integer :: isg,ieg,jsg,jeg  ! global domain


!  ---- version number -----

  character(len=128) :: version = '$Id: fms.F90,v 1.1 2002/01/14 18:22:47 fms Exp $'
  character(len=128) :: tag = '$Name: galway $'

  logical :: do_init = .true.


contains

!#######################################################################
! initializes the fms module/package
! also calls mpp initialization routines and reads fms namelist

subroutine fms_init ( )

 integer :: unit, ierr, io

    if (.not.do_init) return    ! return silently if already called
    do_init = .false.

!---- initialize mpp routines ----

    call mpp_init
    call mpp_domains_init
    call mpp_io_init

!---- read namelist input ----

    call nml_error_init  ! first initialize namelist iostat error codes

    if (file_exist('input.nml')) then
       unit = open_namelist_file ( )
       ierr=1; do while (ierr /= 0)
          read  (unit, nml=fms_nml, iostat=io, end=10)
          ierr = check_nml_error(io,'fms_nml')  ! also initializes nml error codes
       enddo
 10    call mpp_close (unit)
    endif

!---- define mpp stack sizes if non-zero -----

    if (        stack_size > 0) call         mpp_set_stack_size (        stack_size)
    if (domains_stack_size > 0) call mpp_domains_set_stack_size (domains_stack_size)

!---- set severity level for warnings ----

    if ( lowercase(trim(warning_level)) == 'fatal' ) then
            call mpp_set_warn_level ( FATAL )
    else if ( lowercase(trim(warning_level)) == 'warning' ) then
            call mpp_set_warn_level ( WARNING )
    else
            call error_mesg ( 'fms_init',  &
            'invalid entry for namelist variable warning_level', FATAL )
    endif

!--- write version info and namelist to logfile ---

    call write_version_number (version, tag)
    if (mpp_pe() == mpp_root_pe()) then
      write (stdlog(), nml=fms_nml)
      write (stdlog(),*) 'nml_error_codes=', nml_error_codes(1:num_nml_error_codes)
    endif

!---- initialize module domain2d pointer ----

    nullify (Domain)

end subroutine fms_init

!#######################################################################
! terminates the fms module/package
! also calls mpp destructor routines

subroutine fms_end ( )

    if (do_init) return  ! return silently
    call mpp_io_exit
    call mpp_domains_exit
    call mpp_exit
    do_init = .true.

end subroutine fms_end

!#######################################################################
! check the existence of the given file name
! if the file_name string has zero length or the
! first character is blank return a false result

 function file_exist (file_name)
  character(len=*), intent(in) :: file_name
  logical  file_exist

   file_exist = .false.
   if (len_trim(file_name) == 0) return
   if (file_name(1:1) == ' ')    return

   inquire (file=trim(file_name), exist=file_exist)

 end function file_exist

!#######################################################################
! wrapper for the mpp error handler
! users should try to use the mpp_error interface

 subroutine error_mesg (routine, message, level)
  character(len=*), intent(in) :: routine, message
  integer,          intent(in) :: level

!  input:
!      routine   name of the calling routine (character string)
!      message   message written to output   (character string)
!      level     set to NOTE, MESSAGE, or FATAL (integer)

    if (do_init) call fms_init ( )
    call mpp_error ( routine, message, level )

 end subroutine error_mesg

!#######################################################################
! used to check the iostat argument that is
! returned after reading a namelist
! see the online documentation for how this routine might be used

 function check_nml_error (iostat, nml_name) result (error_code)

  integer,          intent(in) :: iostat
  character(len=*), intent(in) :: nml_name
  integer   error_code, i
  character(len=128) :: err_str

   if (do_init) call fms_init ( )

   error_code = iostat

   do i = 1, num_nml_error_codes
        if (error_code == nml_error_codes(i)) return
   enddo

!  ------ fatal namelist error -------
!  ------ only on root pe ----------------
   if (mpp_pe() == mpp_root_pe()) then
       write (err_str,*) 'while reading namelist ',  &
                         trim(nml_name), ', iostat = ',error_code
       call error_mesg ('check_nml_error in fms_mod', err_str, FATAL)
   endif

end function check_nml_error

!-----------------------------------------------------------------------
!   private routine for initializing allowable error codes

subroutine nml_error_init

! some compilers return non-zero iostat values while
! reading through files with multiple namelist records
! this routines "attempts" to identify the iostat values associated
! with records not belonging to the requested namelist

   integer  unit, io, ir
   real    ::  a=1.
   integer ::  b=1
   logical ::  c=.true.
   character(len=8) ::  d='testing'
   namelist /b_nml/  a,b,c,d

      nml_error_codes(1) = 0

!     ---- create dummy namelist file that resembles actual ----
!     ---- (each pe has own copy) ----
      call mpp_open (unit, '_read_error.nml', form=MPP_ASCII,  &
                     action=MPP_OVERWR, access=MPP_SEQUENTIAL, &
                     threading=MPP_MULTI)
!     ---- due to namelist bug this will not always work ---
      write (unit, 10)
  10  format ('    ', &
             /' &a_nml  a=1.  /',    &
             /'#------------------', &
             /' &b_nml  a=5., b=0, c=.false., d=''test'',  &end')
      call mpp_close (unit)

!     ---- read namelist files and save error codes ----
      call mpp_open (unit, '_read_error.nml', form=MPP_ASCII,  &
                     action=MPP_RDONLY, access=MPP_SEQUENTIAL, &
                     threading=MPP_MULTI)
      ir=1; io=1; do
         read  (unit, nml=b_nml, iostat=io, end=20)
         if (io == 0) exit
         ir=ir+1; nml_error_codes(ir)=io
      enddo
  20  call mpp_close (unit, action=MPP_DELETE)

      num_nml_error_codes = ir
!del  if (mpp_pe() == mpp_root_pe()) &
!del  print *, 'PE,nml_error_codes=',mpp_pe(), nml_error_codes(1:ir)
      do_nml_error_init = .false.

end subroutine nml_error_init

!#######################################################################
! prints module version number to the log file of specified unit number

 subroutine write_version_number (version, tagname, unit)

!   in:  version = string that contains routine name and version number
!
!   optional in:
!        tagname = cvs tag name that code was checked out with
!        unit    = alternate unit number to direct output  
!                  (default: unit=stdlog)

   character(len=*), intent(in) :: version
   character(len=*), intent(in), optional :: tagname 
   integer,          intent(in), optional :: unit 

   integer :: logunit 

   if (do_init) call fms_init ( )

     logunit = stdlog()
     if (present(unit)) then
         logunit = unit
     else    
       ! only allow stdlog messages on root pe
         if ( mpp_pe() /= mpp_root_pe() ) return
     endif   

     if (present(tagname)) then
         write (logunit,'(/,80("="),/(a))') trim(version), trim(tagname)
     else    
         write (logunit,'(/,80("="),/(a))') trim(version)
     endif   

 end subroutine write_version_number

!#######################################################################
!#######################################################################
!
! routines for opening specific types of files:
!
!                       form        action 
! open_namelist_file  MPP_ASCII   MPP_RDONLY  
! open restart_file   MPP_NATIVE
! open_ieee32_file    MPP_IEEE32
!
! all have: access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true.
! use the close_file interface to close these files
!
! if other types of files need to be opened the mpp_open and
! mpp_close routines in the mpp_io_mod should be used
!
!#######################################################################
! opens single namelist file for reading only by all PEs
! the default file opened is called "input.nml"

 function open_namelist_file (file) result (unit)
 character(len=*), intent(in), optional :: file
 integer :: unit

   if (do_init) call fms_init ( )

   if (present(file)) then
       call mpp_open ( unit, file, form=MPP_ASCII, action=MPP_RDONLY, &
                       access=MPP_SEQUENTIAL, threading=MPP_SINGLE )
   else
       call mpp_open ( unit, 'input.nml', form=MPP_ASCII, action=MPP_RDONLY, &
                       access=MPP_SEQUENTIAL, threading=MPP_SINGLE )
   endif

 end function open_namelist_file

!#######################################################################
! opens single restart file for reading by all PEs or
! writing by root PE only
! the file has native format and no mpp header records

 function open_restart_file (file, action) result (unit)
 character(len=*), intent(in) :: file, action
 integer :: unit

 integer :: mpp_action

   if (do_init) call fms_init ( )

!   --- action (read,write) ---

    select case (lowercase(trim(action)))
       case ('read')
           mpp_action = MPP_RDONLY
       case ('write')
           mpp_action = MPP_OVERWR
       case default
           call error_mesg ('open_restart_file in fms_mod', &
                            'invalid option for argument action', FATAL)
    end select

    call mpp_open ( unit, file, form=MPP_NATIVE, action=mpp_action, &
            access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true. )

 end function open_restart_file

!#######################################################################
! opens single 32-bit ieee file for reading by all PEs or 
! writing by root PE only (writing is not recommended)
! the file has no mpp header records

 function open_ieee32_file (file, action) result (unit)
 character(len=*), intent(in) :: file, action
 integer :: unit

 integer :: mpp_action

   if (do_init) call fms_init ( )

!   --- action (read,write) ---

    select case (lowercase(trim(action)))
       case ('read')
           mpp_action = MPP_RDONLY
       case ('write')
           mpp_action = MPP_OVERWR
       case default
           call error_mesg ('open_ieee32_file in fms_mod', &
                            'invalid option for argument action', FATAL)
    end select

    if (iospec_ieee32(1:1) == ' ') then
       call mpp_open ( unit, file, form=MPP_IEEE32, action=mpp_action, &
                       access=MPP_SEQUENTIAL, threading=MPP_SINGLE,    &
                       nohdrs=.true. )
    else
       call mpp_open ( unit, file, form=MPP_IEEE32, action=mpp_action, &
                       access=MPP_SEQUENTIAL, threading=MPP_SINGLE,    &
                       nohdrs=.true., iospec=iospec_ieee32 )
    endif

 end function open_ieee32_file

!#######################################################################

 subroutine close_file (unit, status)
   integer,          intent(in)           :: unit
   character(len=*), intent(in), optional :: status

   if (do_init) call fms_init ( )

   if (unit == stdlog()) return

   if (present(status)) then
       if (lowercase(trim(status)) == 'delete') then
           call mpp_close (unit, action=MPP_DELETE)
       else
           call error_mesg ('close_file in fms_mod', &
                            'invalid value for status', FATAL)
       endif
   else
           call mpp_close (unit)
   endif
   
 end subroutine close_file

!#######################################################################
!#######################################################################
! set_domain is called to save the domain2d data type prior to
! calling the distributed data I/O routines, read_data and write_data

 subroutine set_domain (Domain2)

   type(domain2D), intent(in), target :: Domain2

   if (do_init) call fms_init ( )

!  --- set_domain must be called before a read_data or write_data ---

   if (associated(Domain)) nullify (Domain)
   Domain => Domain2

!  --- module indexing to shorten read/write routines ---

   call mpp_get_compute_domain (Domain,is ,ie ,js ,je )
   call mpp_get_data_domain    (Domain,isd,ied,jsd,jed)
   call mpp_get_global_domain  (Domain,isg,ieg,jsg,jeg)

!-----------------------------------------------------------------------

 end subroutine set_domain

!#######################################################################
! this will be a private routine with the next release
! users should get the domain decomposition from the domain2d data type

 subroutine get_domain_decomp ( x, y )

  integer, intent(out), dimension(4) :: x, y

  if (mpp_pe() == mpp_root_pe())  call error_mesg ('fms_mod', &
                        'subroutine get_domain_decomp will be &
                         &removed with the next release', NOTE)
  x = (/ isg, ieg, is, ie /)
  y = (/ jsg, jeg, js, je /)

 end subroutine get_domain_decomp

!#######################################################################
!#######################################################################
!   --------- routines for reading distributed data ---------
! before calling these routines the domain decompostion must be set
! by calling "set_domain" with the appropriate domain2d data type
!
! reading can be done either by all PEs (default) or by only the root PE
! this is controlled by namelist variable "read_all_pe".
!#######################################################################

 subroutine read_data_2d ( unit, data, end )

   integer, intent(in)                        :: unit
   real,    intent(out), dimension(isd:,jsd:) :: data
   logical, intent(out), optional             :: end

   real, dimension(isg:ieg,jsg:jeg) :: gdata
   integer :: len

!------------------------------
include "read_data_2d.inc"
!------------------------------

 end subroutine read_data_2d

!#######################################################################

 subroutine read_ldata_2d ( unit, data, end )

   integer, intent(in)                        :: unit
   logical, intent(out), dimension(isd:,jsd:) :: data
   logical, intent(out), optional             :: end

   logical, dimension(isg:ieg,jsg:jeg) :: gdata
   integer :: len

!------------------------------
include "read_data_2d.inc"
!------------------------------

 end subroutine read_ldata_2d

!#######################################################################

 subroutine read_idata_2d ( unit, data, end )

   integer, intent(in)                        :: unit
   integer, intent(out), dimension(isd:,jsd:) :: data
   logical, intent(out), optional             :: end

   integer, dimension(isg:ieg,jsg:jeg) :: gdata
   integer :: len

!------------------------------
include "read_data_2d.inc"
!------------------------------

 end subroutine read_idata_2d

!#######################################################################

 subroutine read_cdata_2d ( unit, data, end )

   integer, intent(in)                        :: unit
   complex,    intent(out), dimension(isd:,jsd:) :: data
   logical, intent(out), optional             :: end

   complex, dimension(isg:ieg,jsg:jeg) :: gdata
   integer :: len

!------------------------------
 include "read_data_2d.inc"
!------------------------------

 end subroutine read_cdata_2d

!#######################################################################

 subroutine read_data_3d ( unit, data, end )

   integer, intent(in)                          :: unit
   real,    intent(out), dimension(isd:,jsd:,:) :: data
   logical, intent(out), optional               :: end

   real, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata
   integer :: len

!------------------------------
 include "read_data_3d.inc"
!------------------------------

 end subroutine read_data_3d

!#######################################################################

 subroutine read_cdata_3d ( unit, data, end )

   integer, intent(in)                          :: unit
   complex, intent(out), dimension(isd:,jsd:,:) :: data
   logical, intent(out), optional               :: end

   complex, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata
   integer :: len

!------------------------------
 include "read_data_3d.inc"
!------------------------------

 end subroutine read_cdata_3d

!#######################################################################

 subroutine read_data_4d ( unit, data, end )

   integer, intent(in)                            :: unit
   real,    intent(out), dimension(isd:,jsd:,:,:) :: data
   logical, intent(out), optional                 :: end

   real, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
   integer :: len

! WARNING: memory usage with this routine could be costly
!------------------------------
 include "read_data_4d.inc"
!------------------------------

 end subroutine read_data_4d

!#######################################################################

 subroutine read_cdata_4d ( unit, data, end )

   integer, intent(in)                            :: unit
   complex, intent(out), dimension(isd:,jsd:,:,:) :: data
   logical, intent(out), optional                 :: end

   complex, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
   integer :: len

! WARNING: memory usage with this routine could be costly
!------------------------------
 include "read_data_4d.inc"
!------------------------------

 end subroutine read_cdata_4d

!#######################################################################
!     -------- routines for writing distributed data --------
! before calling these routines the domain decompostion must be set
! by calling "set_domain" with the appropriate domain2d data type
!#######################################################################

 subroutine write_data_2d ( unit, data )

   integer, intent(in)                       :: unit
   real,    intent(in), dimension(isd:,jsd:) :: data

   real, dimension(isg:ieg,jsg:jeg) :: gdata

!--------------------------
 include "write_data.inc"
!--------------------------

 end subroutine write_data_2d

!#######################################################################

 subroutine write_ldata_2d ( unit, data )

   integer, intent(in)                       :: unit
   logical, intent(in), dimension(isd:,jsd:) :: data

   logical, dimension(isg:ieg,jsg:jeg) :: gdata

!--------------------------
 include "write_data.inc"
!--------------------------

 end subroutine write_ldata_2d

!#######################################################################

 subroutine write_idata_2d ( unit, data )

   integer, intent(in)                       :: unit
   integer, intent(in), dimension(isd:,jsd:) :: data

   integer, dimension(isg:ieg,jsg:jeg) :: gdata

!--------------------------
 include "write_data.inc"
!--------------------------

end subroutine write_idata_2d

!#######################################################################

 subroutine write_cdata_2d ( unit, data )

   integer, intent(in)                       :: unit
   complex, intent(in), dimension(isd:,jsd:) :: data

   complex, dimension(isg:ieg,jsg:jeg) :: gdata

!--------------------------
 include "write_data.inc"
!--------------------------

 end subroutine write_cdata_2d

!#######################################################################

 subroutine write_data_3d ( unit, data )

   integer, intent(in) :: unit
   real,    intent(in), dimension(isd:,jsd:,:) :: data

   real, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata

!--------------------------
 include "write_data.inc"
!--------------------------

 end subroutine write_data_3d

!#######################################################################

 subroutine write_cdata_3d ( unit, data )

   integer, intent(in) :: unit
   complex, intent(in), dimension(isd:,jsd:,:) :: data

   complex, dimension(isg:ieg,jsg:jeg,size(data,3)) :: gdata

!--------------------------
 include "write_data.inc"
!--------------------------

 end subroutine write_cdata_3d

!#######################################################################

 subroutine write_data_4d ( unit, data )

   integer, intent(in) :: unit
   real,    intent(in), dimension(isd:,jsd:,:,:) :: data

   real, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
   integer :: n

!--------------------------

   if (.not.associated(Domain)) call error_mesg &
        ('write_data in fms_mod', 'set_domain not called', FATAL)

 ! get the global data and write only on root pe
 ! do this one field at a time to save memory
   do n = 1, size(data,4)
    call mpp_global_field ( Domain, data(:,:,:,n), gdata(:,:,:,n) )
   enddo
   if ( mpp_pe() == mpp_root_pe() ) write (unit) gdata

!--------------------------

 end subroutine write_data_4d

!#######################################################################

 subroutine write_cdata_4d ( unit, data )

   integer, intent(in) :: unit
   complex,    intent(in), dimension(isd:,jsd:,:,:) :: data

   complex, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
   integer :: n

!--------------------------

   if (.not.associated(Domain)) call error_mesg &
        ('write_data in fms_mod', 'set_domain not called', FATAL)

 ! get the global data and write only on root pe
 ! do this one field at a time to save memory
   do n = 1, size(data,4)
    call mpp_global_field ( Domain, data(:,:,:,n), gdata(:,:,:,n) )
   enddo
   if ( mpp_pe() == mpp_root_pe() ) write (unit) gdata

!--------------------------

 end subroutine write_cdata_4d

!#######################################################################
! private routines (read_eof,do_read)
! this routine is called when an EOF is found while
! reading a distributed data file using read_data

 subroutine read_eof (end)
   logical, intent(out), optional :: end

   if (present(end))then
       end = .true.
   else
       call error_mesg ('fms_mod', 'unexpected EOF', FATAL)
   endif

 end subroutine read_eof

!#######################################################################
! determines if current pe should read data
! checks namelist variable read_all_pe

 function do_read ( )
 logical :: do_read

       do_read = mpp_pe() == mpp_root_pe() .or. read_all_pe

 end function do_read

!#######################################################################
!#######################################################################
! routines for timing sections of code
!       mpp_clock_init (wrapper for mpp_clock_id)
!       mpp_clock_begin
!       mpp_clock_end
!#######################################################################

 function mpp_clock_init ( name, level, flags ) result (id)
 character(len=*),  intent(in) :: name
 integer,           intent(in) :: level
 integer, optional, intent(in) :: flags
 integer                       :: id

    if (do_init) call fms_init ( )

  ! only register this clock when "timing_level"
  ! is .GE. then this clock's (timing) level
  ! otherwise return a zero id

    if ( level <= timing_level ) then
        id = mpp_clock_id (name, flags)
    else
        id = 0
    endif

 end function mpp_clock_init

!#######################################################################
!  functions for changing the case of character strings
!#######################################################################

!   change to all lower case

 function lowercase (cs) 
 character(len=*), intent(in) :: cs
 character(len=len(cs))       :: lowercase 
 character :: ca(len(cs)) 

 integer, parameter :: co=iachar('a')-iachar('A') ! case offset
    
    ca = transfer(cs,"x",len(cs)) 
    where (ca >= "A" .and. ca <= "Z") ca = achar(iachar(ca)+co) 
    lowercase = transfer(ca,cs) 
    
 end function lowercase 

!#######################################################################

!   change to all upper case

 function uppercase (cs) 
 character(len=*), intent(in) :: cs
 character(len=len(cs))       :: uppercase 
 character :: ca(len(cs)) 

 integer, parameter :: co=iachar('A')-iachar('a') ! case offset
    
    ca = transfer(cs,"x",len(cs)) 
    where (ca >= "a" .and. ca <= "z") ca = achar(iachar(ca)+co) 
    uppercase = transfer(ca,cs) 
    
 end function uppercase 

!#######################################################################

! match the input character string to a string
! in an array/list of character strings

function string_array_index ( string, string_array, index ) result (found)
character(len=*),  intent(in)  :: string, string_array(:)
integer, optional, intent(out) :: index
logical :: found
integer :: i

! initialize this function to false
! loop thru string_array and exit when a match is found

  found = .false.
  if (present(index)) index = 0

  do i = 1, size(string_array)
    ! found a string match ?
    if ( trim(string) == trim(string_array(i)) ) then
         found = .true.
         if (present(index)) index = i
         exit
    endif
  enddo

end function string_array_index

!#######################################################################

! determines if the real input array has
! monotonically increasing or decreasing values

function monotonic_array ( array, direction )
real,    intent(in)            :: array(:)
integer, intent(out), optional :: direction
logical :: monotonic_array
integer :: i

! initialize
  monotonic_array = .false.
  if (present(direction)) direction = 0

! array too short
  if ( size(array) < 2 ) return

! ascending
  if ( array(1) < array(size(array)) ) then
     do i = 2, size(array)
       if (array(i-1) < array(i)) cycle
       return
     enddo
     monotonic_array = .true.
     if (present(direction)) direction = +1

! descending
  else
     do i = 2, size(array)
       if (array(i-1) > array(i)) cycle
       return
     enddo
     monotonic_array = .true.
     if (present(direction)) direction = -1
  endif

end function monotonic_array

!#######################################################################

end module fms_mod

