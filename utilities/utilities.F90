
module utilities_mod

!-----------------------------------------------------------------------
!
!   A collection of simple useful programs.
!
!    file_exist       Function that returns if a given
!                     file name exists
!
!    error_mesg       Print warning and error messages, 
!                     terminates program for error messages.
!                     (use error levels NOTE,WARNING,FATAL)
!
!    check_nml_error  Checks the iostat returned when reading
!                     namelists and determines if the error code
!                     is valid, if not the program is terminated.
!
!    open_file        Opens the given file name for i/o and returns
!                     a unit number.  If the file is already open
!                     the unit number is returned.
!
!    close_file       Closes the given unit number for i/o.
!
!    print_version_number    Prints out a routine name and
!                            version number to a specified unit
!
!    set_domain          Sets a pointer for a private utilities_mod
!                        domain2d data type. This domain2d type is
!                        used for subsequent calls by read_data and
!                        write data.
!
!    get_domain_decomp   Returns domain information from the domain2d
!                        type set by the most recent call to set_domain.
!
!    read_data           Reads distributed data.
!
!    write_data          Writes distributed data.
!
!    set_system_clock    Sets (or resets) a real-time clock.
!                        This is automatically called by utilities_init.
!
!    check_system_clock  Outputs the elapsed real-time since the last
!                        call to set_system_clock (or utilities_init).
!
!    utilities_init      Initializes the utilities module and also the
!                        mpp_io module (which initializes all mpp mods).
!                        Will be called automatically if the user does
!                        not call it.
!
!    utilities_end       Calls mpp_io_exit and check_system_clock.
!
!    lowercase           Convert character strings to all lower case
!
!    uppercase           Convert character strings to all upper case
!
!    mpp_clock_init      Sets up a identifier for performance timing
!                          (similar to mpp_clock_id)
!
!-----------------------------------------------------------------------
!  ---- Routines with exact functionality to mpp_mod routines ----
!
!    get_my_pe           Returns the pe number. (mpp_pe)
!
!    get_num_pes         Returns the total number pe's. (mpp_npes)
!
!    sync_all_pes        Synchronizes pe's. (mpp_sync)
!
!    mpp_sum             Compute summation. (mpp_sum)
!    mpp_min             Find minimum. (mpp_min)
!    mpp_max             Find maximum. (mpp_max)
!
!    check_sum           Performs integer check sums (mpp_chksum)
!
!    mpp_clock_begin     Initiates performance timing (mpp_clock_begin)
!    mpp_clock_end       Terminates performance timing (mpp_clock_end)
!
!-----------------------------------------------------------------------

use          mpp_mod, only:  mpp_error, NOTE, WARNING, FATAL, &
                             mpp_set_warn_level,              &
                             mpp_transmit, ALL_PES,           &
                             mpp_sum, mpp_min, mpp_max,       &
                             get_my_pe => mpp_pe,             &
                             get_num_pes => mpp_npes,         &
                             sync_all_pes => mpp_sync,        &
                             check_sum => mpp_chksum,         &
                             mpp_clock_begin, mpp_clock_end,  &
                             mpp_clock_id, mpp_exit,          &
                             mpp_set_stack_size

use  mpp_domains_mod, only:  domain2D, mpp_define_domains, &
                             mpp_update_domains, GLOBAL_DATA_DOMAIN, &
                             mpp_get_global, mpp_domains_exit,       &
                             mpp_domains_set_stack_size

use       mpp_io_mod, only:  mpp_io_init, mpp_open, mpp_close,         &
                       MPP_ASCII, MPP_NATIVE, MPP_IEEE32, MPP_NETCDF,  &
                       MPP_RDONLY, MPP_WRONLY, MPP_APPEND, MPP_OVERWR, &
                       MPP_SEQUENTIAL, MPP_DIRECT,                     &
                       MPP_SINGLE, MPP_MULTI, MPP_DELETE, mpp_io_exit

implicit none
private

public :: file_exist, open_file, read_data, write_data,  &
          set_domain, get_domain_decomp,                 &
          check_nml_error, print_version_number,         &
          error_mesg, NOTE, WARNING, FATAL,              &
          get_my_pe, get_num_pes, sync_all_pes,          &
          mpp_sum, mpp_min, mpp_max,                     &
          utilities_init, utilities_end, close_file,     &
          set_system_clock, check_system_clock, check_sum, &
          mpp_clock_begin, mpp_clock_end, mpp_clock_init,  &
          lowercase, uppercase

!   ---- private data for check_nml_error ----

   integer, private :: num_nml_error_codes, nml_error_codes(20)
   logical, private :: do_nml_error_init = .true.
   private  nml_error_init


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

  logical           :: time_all_pe   = .false.
  integer           :: timing_level  = 0
  logical           :: read_all_pe   = .true.
  character(len=8)  :: warning_level = 'fatal'
  character(len=64) :: iospec_ieee32 = '-F f77,cachea:48:1'
  logical           :: one_level_restarts = .false.
  integer           :: stack_size = 0
  integer           :: domains_stack_size = 0

  namelist /utilities_nml/  timing_level, read_all_pe, warning_level, &
                            iospec_ieee32, one_level_restarts,        &
                            stack_size, domains_stack_size,           &
                            time_all_pe


!------ private data, pointer to current 2d domain ------

  type(domain2D), pointer :: Domain

  integer :: is,ie,js,je, isd,ied,jsd,jed, isg,ieg,jsg,jeg

!------ private data for system clock routines ------

  integer :: count_init, count_rate, count_max


  integer :: mxdim3  ! maximum number of levels written/read per record

!  ---- version number -----

  character(len=128) :: version = '$Id: utilities.F90,v 1.4 2001/03/06 21:05:30 fms Exp $'
  character(len=128) :: tag = '$Name: damascus $'

  logical :: do_init = .true.


contains

!#######################################################################

   function file_exist (file_name)

      character(len=*), intent(in) :: file_name
      logical  file_exist

      inquire (file=trim(file_name), exist=file_exist)

   end function file_exist

!#######################################################################

   subroutine error_mesg (routine, message, level)

!             ------------------------------------
!             |                                  |
!             |    a very simple error handler   |
!             |                                  |
!             ------------------------------------
!
!  input:
!      routine   name of the calling routine (character string)
!      message   message written to standard output (character string)
!      level     set to NOTE, MESSAGE, or FATAL
!
      character(len=*), intent(in) :: routine, message
      integer,          intent(in) :: level

      character(len=2) :: sep = ', '

      if ( do_init ) call utilities_init ( )

      call mpp_error ( level, trim(routine) // sep // trim(message) )

!-----------------------------------------------------------------------

   end subroutine error_mesg

!#######################################################################

function check_nml_error (iostat, nml_name) result (error_code)

   integer,          intent(in) :: iostat
   character(len=*), intent(in) :: nml_name
   integer   error_code, i
   character(len=128) :: err_str

   if (do_nml_error_init) call nml_error_init

   error_code = iostat

   do i = 1, num_nml_error_codes
        if (error_code == nml_error_codes(i)) return
   enddo

!  ------ fatal namelist error -------
!  ------ only on pe0 ----------------
   if (get_my_pe() == 0) then
       write (err_str,*) 'while reading namelist ',  &
                         trim(nml_name), ', iostat = ',error_code
       call error_mesg ('check_nml_error', err_str, FATAL)
   endif

end function check_nml_error

!-----------------------------------------------------------------------

subroutine nml_error_init

!   private routine for initializing allowable error codes

   integer  unit, io, ir
   real    ::  a=1.
   integer ::  b=1
   logical ::  c=.true.
   character(len=8) ::  d='testing'
   namelist /b_nml/  a,b,c,d

      nml_error_codes(1) = 0

!     ---- create dummy namelist file that resembles actual ----
!     ---- (each pe has own copy) ----
      unit = open_file ('_read_error.nml', action='write',  &
                                           threading='multi')
!     ---- due to namelist bug this will not always work ---
      write (unit, 10)
  10  format ('    ', &
             /' &a_nml  a=1.  /',    &
             /'#------------------', &
             /' &b_nml  e=5.,  &end')
!bug         /' &b_nml  a=5., b=0, c=.false., d=''test'',  &end')
      call close_file (unit)

!     ---- read namelist files and save error codes ----
      unit = open_file ('_read_error.nml', action='read',  &
                                           threading='multi')
      ir=1; io=1; do
         read  (unit, nml=b_nml, iostat=io, end=20)
         if (io == 0) exit
         ir=ir+1; nml_error_codes(ir)=io
      enddo
  20  call close_file (unit, status='delete')

      num_nml_error_codes = ir
!del  if (get_my_pe() == 0) &
!del  print *, 'PE,nml_error_codes=',get_my_pe(), nml_error_codes(1:ir)
      do_nml_error_init = .false.

end subroutine nml_error_init

!#######################################################################
!#######################################################################

   function open_file (file, form, action, access, threading, recl) &
               result (unit)

   character(len=*), intent(in) :: file
   character(len=*), intent(in), optional :: form, action, access, &
                                             threading
   integer         , intent(in), optional :: recl
   integer  :: unit

   integer           :: nc
   character(len=11) :: form11
   character(len=10) :: access10
   character(len=6)  :: action6, thread6
   character(len=64) :: iospec
   integer :: mpp_format, mpp_action, mpp_access, mpp_thread
   logical :: open, no_headers
!-----------------------------------------------------------------------

   if ( do_init ) call utilities_init ( )

!   ---- is this file open and connected to a unit ?? ----

   inquire (file=trim(file), opened=open, number=unit)

   if ( open .and. unit >= 0 ) return

!   -------------- default -----------


!   -------------- file format -------------------

         mpp_format = MPP_ASCII
         iospec = ' '

         if (present(form)) then
             nc = min(11,len(form))
             form11(1:nc) = form(1:nc)

             if ( form11(1:9) == 'formatted' .or. &
                  form11(1:9) == 'FORMATTED' .or. &
                  form11(1:5) == 'ascii'     .or. &
                  form11(1:5) == 'ASCII' )   then
                    mpp_format = MPP_ASCII

             else if ( form11(1:11) == 'unformatted' .or. &
                       form11(1:11) == 'UNFORMATTED' .or. &
                       form11(1:6)  == 'native'      .or. &
                       form11(1:6)  == 'NATIVE' )    then
                    mpp_format = MPP_NATIVE

             else if ( form11(1:6) == 'ieee32' .or. &
                       form11(1:6) == 'IEEE32' ) then
                    mpp_format = MPP_IEEE32
                    iospec     = iospec_ieee32

             else if ( form11(1:6) == 'netcdf' .or. &
                       form11(1:6) == 'NETCDF' ) then
                    mpp_format = MPP_NETCDF

             endif

         endif

!   ----------------------------------------------

         if ( present(action) ) then
              nc = min(6,len(action))
              action6(1:nc) = action(1:nc)

              if ( action6(1:4) == 'read' .or. &
                   action6(1:4) == 'READ' ) then
                      mpp_action = MPP_RDONLY

              else if ( action6(1:5) == 'write' .or. &
                        action6(1:5) == 'WRITE' ) then
                      mpp_action = MPP_OVERWR

              else if ( action6(1:6) == 'append' .or. &
                        action6(1:6) == 'APPEND' ) then
                      mpp_action = MPP_APPEND

              endif

         endif
             
!   ----------------------------------------------

         mpp_access = MPP_SEQUENTIAL

         if ( present(access) ) then
              nc = min(10,len(access))
              access10(1:nc) = access(1:nc)

              if ( access10(1:6) == 'direct' .or. &
                   access10(1:6) == 'DIRECT' ) then
                      mpp_access = MPP_DIRECT

              else if ( access10(1:10) == 'sequential' .or. &
                        access10(1:10) == 'SEQUENTIAL' ) then
                      mpp_access = MPP_SEQUENTIAL

              endif

         endif
        
!   ----------------------------------------------

         mpp_thread = MPP_SINGLE
         no_headers = .true.

         if ( present(threading) ) then
              nc = min(6,len(threading))
              thread6(1:nc) = threading(1:nc)

              if ( thread6(1:6) == 'single' .or. &
                   thread6(1:6) == 'SINGLE' ) then
                      mpp_thread = MPP_SINGLE
                      no_headers = .true.

              else if ( thread6(1:5) == 'multi' .or. &
                        thread6(1:5) == 'MULTI' ) then
                      mpp_thread = MPP_MULTI

                      if (trim(file) /= '_read_error.nml') &
                      no_headers = .false.

              endif

         endif
        
!   ----------------------------------------------

 if ( iospec(1:1) == ' ' ) then
    call mpp_open ( unit, file, form=mpp_format, action=mpp_action, &
                    access=mpp_access, threading=mpp_thread,        &
                                   nohdrs=no_headers, recl=recl )
 else
    call mpp_open ( unit, file, form=mpp_format, action=mpp_action, &
                    access=mpp_access, threading=mpp_thread,        &
                    iospec=iospec, nohdrs=no_headers, recl=recl )
 endif

!-----------------------------------------------------------------------

   end function open_file

!#######################################################################

   subroutine print_version_number (unit, routine, version)

! *** prints routine name and version number to a log file ***
!
!    in:  unit    = unit number to direct output
!         routine = routine name (character, max len = 20)
!         version = version name or number (character, max len = 8)

   integer,          intent(in) :: unit
   character(len=*), intent(in) :: routine, version

   integer           :: n
   character(len=20) :: name
   character(len=8)  :: vers

     if ( get_my_pe() /= 0 ) return

     n = min(len(routine),20); name = adjustl(routine(1:n))
     n = min(len(version), 8); vers = adjustl(version(1:n))

     if (unit > 0) then
         write (unit,10) name, vers
     else
         write (*,10) name, vers
     endif

  10 format (/,66('-'),  &
             /,10x, 'ROUTINE = ',a20, '  VERSION = ', a8, &
             /,66('-'))

! 10 format (/,1x, 12('>'), 1x, 'ROUTINE = ',a20, '  VERSION = ', a8, &
!              1x, 12('<'),/)

   end subroutine print_version_number

!#######################################################################

subroutine set_domain (Domain2)

   type(domain2D), intent(in), target :: Domain2

!  --- set_domain must be called before a read_data or write_data ---

   if (associated(Domain)) nullify (Domain)
   Domain => Domain2

!  --- module indexing to shorten read/write routines ---

   is = Domain%X%Compute%start_index
   ie = Domain%X%Compute%end_index
   js = Domain%Y%Compute%start_index
   je = Domain%Y%Compute%end_index

   isd = Domain%X%Data%start_index
   ied = Domain%X%Data%end_index
   jsd = Domain%Y%Data%start_index
   jed = Domain%Y%Data%end_index

   isg = Domain%X%Global%start_index
   ieg = Domain%X%Global%end_index
   jsg = Domain%Y%Global%start_index
   jeg = Domain%Y%Global%end_index

!-----------------------------------------------------------------------

end subroutine set_domain

!#######################################################################

subroutine get_domain_decomp ( x, y )

  integer, intent(out), dimension(4) :: x, y

  x = (/ isg, ieg, is, ie /)
  y = (/ jsg, jeg, js, je /)

end subroutine get_domain_decomp

!#######################################################################

subroutine read_data_2d ( unit, data, end )

   integer, intent(in)                        :: unit
   real,    intent(out), dimension(isd:,jsd:) :: data
   logical, intent(out), optional             :: end

   real, dimension(isg:ieg,jsg:jeg) :: gdata
   logical :: do_read
   integer :: len

!-----------------------------------------------------------------------

   if (.not.associated(Domain)) call error_mesg &
        ('read_data in utilities_mod', 'set_domain not called', FATAL)

       if (present(end)) end = .false.

       do_read = get_my_pe() == 0 .or. read_all_pe

       if (do_read) read (unit,end=10) gdata
       if (.not.read_all_pe) then
           len = size(gdata,1)*size(gdata,2)
           call mpp_transmit ( gdata(isg,jsg), len, ALL_PES, &
                               gdata(isg,jsg), len, 0        )
       endif
       data(is:ie,js:je) = gdata(is:ie,js:je)
       return

   10  if (present(end))then
           end = .true.
       else
           call error_mesg ('utilities_mod', 'unexpected EOF', FATAL)
       endif

!-----------------------------------------------------------------------

end subroutine read_data_2d

!#######################################################################

subroutine read_ldata_2d ( unit, data, end )

   integer, intent(in)                        :: unit
   logical, intent(out), dimension(isd:,jsd:) :: data
   logical, intent(out), optional             :: end

   logical, dimension(isg:ieg,jsg:jeg) :: gdata

!-----------------------------------------------------------------------

   if (.not.associated(Domain)) call error_mesg &
        ('read_data in utilities_mod', 'set_domain not called', FATAL)

       if (present(end)) end = .false.

       read (unit,end=10) gdata
       data(is:ie,js:je) = gdata(is:ie,js:je)
       return

   10  if (present(end))then
           end = .true.
       else
           call error_mesg ('utilities_mod', 'unexpected EOF', FATAL)
       endif

!-----------------------------------------------------------------------

end subroutine read_ldata_2d

!#######################################################################

subroutine read_idata_2d ( unit, data, end )

   integer, intent(in)                        :: unit
   integer, intent(out), dimension(isd:,jsd:) :: data
   logical, intent(out), optional             :: end

   integer, dimension(isg:ieg,jsg:jeg) :: gdata
   logical :: do_read
   integer :: len

!-----------------------------------------------------------------------

   if (.not.associated(Domain)) call error_mesg &
        ('read_data in utilities_mod', 'set_domain not called', FATAL)

       if (present(end)) end = .false.

       do_read = get_my_pe() == 0 .or. read_all_pe

       if (do_read) read (unit,end=10) gdata
       if (.not.read_all_pe) then
           len = size(gdata,1)*size(gdata,2)
           call mpp_transmit ( gdata(isg,jsg), len, ALL_PES, &
                               gdata(isg,jsg), len, 0        )
       endif
       data(is:ie,js:je) = gdata(is:ie,js:je)
       return

   10  if (present(end))then
           end = .true.
       else
           call error_mesg ('utilities_mod', 'unexpected EOF', FATAL)
       endif

!-----------------------------------------------------------------------

end subroutine read_idata_2d

!#######################################################################

subroutine read_cdata_2d ( unit, data, end )

   integer, intent(in)                        :: unit
   complex,    intent(out), dimension(isd:,jsd:) :: data
   logical, intent(out), optional             :: end

   complex, dimension(isg:ieg,jsg:jeg) :: gdata
   logical :: do_read
   integer :: len

!-----------------------------------------------------------------------

   if (.not.associated(Domain)) call error_mesg &
        ('read_data in utilities_mod', 'set_domain not called', FATAL)

       if (present(end)) end = .false.

       do_read = get_my_pe() == 0 .or. read_all_pe

       if (do_read) read (unit,end=10) gdata
       if (.not.read_all_pe) then
           len = size(gdata,1)*size(gdata,2)
           call mpp_transmit ( gdata(isg,jsg), len, ALL_PES, &
                               gdata(isg,jsg), len, 0        )
       endif
       data(is:ie,js:je) = gdata(is:ie,js:je)
       return

   10  if (present(end))then
           end = .true.
       else
           call error_mesg ('utilities_mod', 'unexpected EOF', FATAL)
       endif

!-----------------------------------------------------------------------

end subroutine read_cdata_2d

!#######################################################################

subroutine read_data_3d ( unit, data, end )

   integer, intent(in)                          :: unit
   real,    intent(out), dimension(isd:,jsd:,:) :: data
   logical, intent(out), optional               :: end

   real, dimension(isg:ieg,jsg:jeg,min(size(data,3), mxdim3)) :: gdata
   integer :: m, m1, m2, msize
   logical :: do_read
   integer :: len

!-----------------------------------------------------------------------

   if (.not.associated(Domain)) call error_mesg &
        ('read_data in utilities_mod', 'set_domain not called', FATAL)

       if (present(end)) end = .false.

       do_read = get_my_pe() == 0 .or. read_all_pe

       msize = min(size(data,3), mxdim3)
       m2 = 0

       do m=1,size(data,3)/msize
	 m1 = m2 + 1
	 m2 = m2 + msize
       if (do_read) read (unit,end=10) gdata
       if (.not.read_all_pe) then
           len = size(gdata,1)*size(gdata,2)*size(gdata,3)
           call mpp_transmit ( gdata(isg,jsg,1), len, ALL_PES, &
                               gdata(isg,jsg,1), len, 0        )
       endif
       data(is:ie,js:je,m1:m2) = gdata(is:ie,js:je,:)

       end do
       return

   10  if (present(end))then
           end = .true.
       else
           call error_mesg ('utilities_mod', 'unexpected EOF', FATAL)
       endif

!-----------------------------------------------------------------------

end subroutine read_data_3d

!#######################################################################

subroutine read_cdata_3d ( unit, data, end )

   integer, intent(in)                          :: unit
   complex,    intent(out), dimension(isd:,jsd:,:) :: data
   logical, intent(out), optional               :: end

   complex, dimension(isg:ieg,jsg:jeg,min(size(data,3),mxdim3)) :: gdata
   integer :: m, m1, m2, msize
   logical :: do_read
   integer :: len

!-----------------------------------------------------------------------

   if (.not.associated(Domain)) call error_mesg &
        ('read_data in utilities_mod', 'set_domain not called', FATAL)

       if (present(end)) end = .false.

       do_read = get_my_pe() == 0 .or. read_all_pe

       msize = min(size(data,3), mxdim3)
       m2 = 0

       do m=1,size(data,3)/msize
	 m1 = m2 + 1
	 m2 = m2 + msize
       if (do_read) read (unit,end=10) gdata
       if (.not.read_all_pe) then
           len = size(gdata,1)*size(gdata,2)*size(gdata,3)
           call mpp_transmit ( gdata(isg,jsg,1), len, ALL_PES, &
                               gdata(isg,jsg,1), len, 0        )
       endif
       data(is:ie,js:je,m1:m2) = gdata(is:ie,js:je,:)

       end do
       return

   10  if (present(end))then
           end = .true.
       else
           call error_mesg ('utilities_mod', 'unexpected EOF', FATAL)
       endif

!-----------------------------------------------------------------------

end subroutine read_cdata_3d

!#######################################################################

subroutine read_data_4d ( unit, data, end )

   integer, intent(in)                            :: unit
   real,    intent(out), dimension(isd:,jsd:,:,:) :: data
   logical, intent(out), optional                 :: end

   real, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
   logical :: do_read
   integer :: len

!-----------------------------------------------------------------------
!**** WARNING: memory usage with this routine could be costly *****

   if (.not.associated(Domain)) call error_mesg &
        ('read_data in utilities_mod', 'set_domain not called', FATAL)

       if (present(end)) end = .false.

       do_read = get_my_pe() == 0 .or. read_all_pe

       if (do_read) read (unit,end=10) gdata
       if (.not.read_all_pe) then
           len = size(gdata,1)*size(gdata,2)*size(gdata,3)*size(gdata,4)
           call mpp_transmit ( gdata(isg,jsg,1,1), len, ALL_PES, &
                               gdata(isg,jsg,1,1), len, 0        )
       endif
       data(is:ie,js:je,:,:) = gdata(is:ie,js:je,:,:)
       return

   10  if (present(end))then
           end = .true.
       else
           call error_mesg ('utilities_mod', 'unexpected EOF', FATAL)
       endif

!-----------------------------------------------------------------------

end subroutine read_data_4d

!#######################################################################

subroutine read_cdata_4d ( unit, data, end )

   integer, intent(in)                            :: unit
   complex,    intent(out), dimension(isd:,jsd:,:,:) :: data
   logical, intent(out), optional                 :: end

   complex, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
   logical :: do_read
   integer :: len

!-----------------------------------------------------------------------
!**** WARNING: memory usage with this routine could be costly *****

   if (.not.associated(Domain)) call error_mesg &
        ('read_data in utilities_mod', 'set_domain not called', FATAL)

       if (present(end)) end = .false.

       do_read = get_my_pe() == 0 .or. read_all_pe

       if (do_read) read (unit,end=10) gdata
       if (.not.read_all_pe) then
           len = size(gdata,1)*size(gdata,2)*size(gdata,3)*size(gdata,4)
           call mpp_transmit ( gdata(isg,jsg,1,1), len, ALL_PES, &
                               gdata(isg,jsg,1,1), len, 0        )
       endif
       data(is:ie,js:je,:,:) = gdata(is:ie,js:je,:,:)
       return

   10  if (present(end))then
           end = .true.
       else
           call error_mesg ('utilities_mod', 'unexpected EOF', FATAL)
       endif

!-----------------------------------------------------------------------

end subroutine read_cdata_4d

!#######################################################################

subroutine write_data_2d ( unit, data )

   integer, intent(in)                       :: unit
   real,    intent(in), dimension(isd:,jsd:) :: data

   real, dimension(isg:ieg,jsg:jeg) :: gdata

   if (.not.associated(Domain)) call error_mesg &
        ('write_data in utilities_mod', 'set_domain not called', FATAL)

!---- put field onto global domain ----

   call mpp_get_global ( Domain, data, gdata )

   if ( get_my_pe() == 0 ) write (unit) gdata

!-----------------------------------------------------------------------

end subroutine write_data_2d

!#######################################################################

subroutine write_ldata_2d ( unit, data )

   integer, intent(in)                       :: unit
   logical, intent(in), dimension(isd:,jsd:) :: data

   real, dimension(isg:ieg,jsg:jeg) :: gdata

   real,    dimension(isd:ied,jsd:jed) :: rdata
   logical, dimension(isg:ieg,jsg:jeg) :: ldata

   if (.not.associated(Domain)) call error_mesg &
        ('write_data in utilities_mod', 'set_domain not called', FATAL)

   where (data)
     rdata = 1.
   elsewhere
     rdata = 0.
   endwhere

!---- put field onto global domain ----

   call mpp_get_global ( Domain, rdata, gdata )

!---- switch back to logical ----

   ldata = ( gdata > 0.99 )

   if ( get_my_pe() == 0 ) write (unit) ldata

!-----------------------------------------------------------------------

end subroutine write_ldata_2d

!#######################################################################

subroutine write_idata_2d ( unit, data )

   integer, intent(in)                       :: unit
   integer, intent(in), dimension(isd:,jsd:) :: data

   integer, dimension(isg:ieg,jsg:jeg) :: gdata

   if (.not.associated(Domain)) call error_mesg &
        ('write_data in utilities_mod', 'set_domain not called', FATAL)

!---- put field onto global domain ----

   call mpp_get_global ( Domain, data, gdata )

   if ( get_my_pe() == 0 ) write (unit) gdata

!-----------------------------------------------------------------------

end subroutine write_idata_2d

!#######################################################################

subroutine write_cdata_2d ( unit, data )

   integer, intent(in)                       :: unit
   complex,    intent(in), dimension(isd:,jsd:) :: data

   complex, dimension(isg:ieg,jsg:jeg) :: gdata

   if (.not.associated(Domain)) call error_mesg &
        ('write_data in utilities_mod', 'set_domain not called', FATAL)

!---- put field onto global domain ----

   call mpp_get_global ( Domain, data, gdata )

   if ( get_my_pe() == 0 ) write (unit) gdata

!-----------------------------------------------------------------------

end subroutine write_cdata_2d

!#######################################################################

subroutine write_data_3d ( unit, data )

   integer, intent(in) :: unit
   real,    intent(in), dimension(isd:,jsd:,:) :: data

   real, dimension(isg:ieg,jsg:jeg, min(size(data,3), mxdim3)) :: gdata
   integer   ::  m, m1, m2, msize

   if (.not.associated(Domain)) call error_mesg &
        ('write_data in utilities_mod', 'set_domain not called', FATAL)

   msize = min (size(data,3), mxdim3)
   m2 = 0

   do m = 1, size(data,3)/msize
     m1 = m2 + 1
     m2 = m2 + msize

!---- put field onto global domain ----

     call mpp_get_global ( Domain, data(:,:,m1:m2), gdata )

     if ( get_my_pe() == 0 ) write (unit) gdata
   end do

!-----------------------------------------------------------------------

end subroutine write_data_3d

!#######################################################################

subroutine write_cdata_3d ( unit, data )

   integer, intent(in) :: unit
   complex,    intent(in), dimension(isd:,jsd:,:) :: data

   complex, dimension(isg:ieg,jsg:jeg,min(size(data,3), mxdim3)) :: gdata
   integer   ::  m, m1, m2, msize

   if (.not.associated(Domain)) call error_mesg &
        ('write_data in utilities_mod', 'set_domain not called', FATAL)

   msize = min (size(data,3), mxdim3)
   m2 = 0

   do m = 1, size(data,3)/msize
     m1 = m2 + 1
     m2 = m2 + msize
!---- put field onto global domain ----

     call mpp_get_global ( Domain, data(:,:,m1:m2), gdata )

     if ( get_my_pe() == 0 ) write (unit) gdata
   end do

!-----------------------------------------------------------------------

end subroutine write_cdata_3d

!#######################################################################

subroutine write_data_4d ( unit, data )

   integer, intent(in) :: unit
   real,    intent(in), dimension(isd:,jsd:,:,:) :: data

   real, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
   integer :: n

   if (.not.associated(Domain)) call error_mesg &
        ('write_data in utilities_mod', 'set_domain not called', FATAL)

!---- put field onto global domain ----

   do n = 1, size(data,4)
    call mpp_get_global ( Domain, data(:,:,:,n), gdata(:,:,:,n) )
   enddo

   if ( get_my_pe() == 0 ) write (unit) gdata

!-----------------------------------------------------------------------

end subroutine write_data_4d

!#######################################################################

subroutine write_cdata_4d ( unit, data )

   integer, intent(in) :: unit
   complex,    intent(in), dimension(isd:,jsd:,:,:) :: data

   complex, dimension(isg:ieg,jsg:jeg,size(data,3),size(data,4)) :: gdata
   integer :: n

   if (.not.associated(Domain)) call error_mesg &
        ('write_data in utilities_mod', 'set_domain not called', FATAL)

!---- put field onto global domain ----

   do n = 1, size(data,4)
    call mpp_get_global ( Domain, data(:,:,:,n), gdata(:,:,:,n) )
   enddo

   if ( get_my_pe() == 0 ) write (unit) gdata

!-----------------------------------------------------------------------

end subroutine write_cdata_4d

!#######################################################################

subroutine set_system_clock

 call system_clock (count_init, count_rate, count_max)

end subroutine set_system_clock

!#######################################################################

subroutine check_system_clock (string)

character(len=*), intent(in), optional :: string

!-----------------------------------------------------------------------
! prints out the elapsed time since the last call to check_system_clock
!-----------------------------------------------------------------------

character(len=51) :: fmt
integer :: nd, nc, np, count, clocks
real    :: seconds

 if ( .not. time_all_pe .and. get_my_pe() /= 0 ) return

!if ( get_my_pe() == 0 ) then
!     call error_mesg ('utilities_mod',  &
!                      'routine check_system_clock is obsolete', NOTE)
!endif

 call system_clock (count)

 clocks  = count - count_init
 seconds = real(clocks)/real(count_rate)

! --- formatting ---

   nd=10-int(log10(seconds));            nd=min(nd,9)
   nc=int(log10(real(clocks)))+1;        nc=max(nc,10)
   np=int(log10(real(get_my_pe()+1)))+1; np=min(np,9)


! --- output to stdout ---

   if (present(string)) then
      write (fmt,10) np,nd,nc
      write (*,fmt) trim(string), get_my_pe(), seconds, clocks
   else
      write (fmt(1:47),11) np,nd,nc
      write (*,fmt(1:47)) get_my_pe(), seconds, clocks
   endif

! --- reset clock ---

   call system_clock (count_init)

10 format('(a,'': pe='',i',i1,',4x,''seconds='',f12.',i1,',4x,''clocks='',i',i2,')')
11 format('(''pe='',i',i1,',4x,''seconds='',f12.',i1,',4x,''clocks='',i',i2,')')
!-----------------------------------------------------------------------

end subroutine check_system_clock

!#######################################################################

 function mpp_clock_init ( module, routine, level ) result (id)
 character(len=*), intent(in) :: module, routine
 integer,          intent(in) :: level
 integer                      :: id
 character(len=128) :: label

    if ( do_init ) call utilities_init ( )

    if ( level == timing_level ) then
        label = trim(routine) // '_in_' // uppercase(trim(module))
        id = mpp_clock_id (trim(label(1:24)))
    else
        id = 0
    endif

 end function mpp_clock_init

!#######################################################################

subroutine utilities_init

  integer :: unit, ierr, io

     if ( .not. do_init ) return
     do_init = .false.

     call set_system_clock

!---- initialize mpp routines ----

     call mpp_io_init

!---- read namelist input ----

    if (file_exist('input.nml')) then
       unit = open_file ('input.nml', action='read')
       ierr=1; do while (ierr /= 0)
          read  (unit, nml=utilities_nml, iostat=io, end=10)
          ierr = check_nml_error(io,'utilities_nml')
       enddo
 10    call close_file (unit)
    endif

!---- define mpp stack sizes if non-zero -----

    if (        stack_size > 0) call         mpp_set_stack_size (        stack_size)
    if (domains_stack_size > 0) call mpp_domains_set_stack_size (domains_stack_size)

!---- define size of third dimension in arrays used when 
!---- reading/writing restarts

    if (one_level_restarts) then
       mxdim3 = 1
    else 
       mxdim3 = 1000000       ! arbitrary, larger than any real dimen
    endif


!---- initialize namelist error codes if necessary ----

    if (do_nml_error_init) call nml_error_init

!---- initialize module domain2d pointer ----

    nullify (Domain)

!---- set severity level for warnings ----

    if ( warning_level(1:5) == 'FATAL' .or.  &
         warning_level(1:5) == 'fatal' ) then
            call mpp_set_warn_level ( FATAL )
    else if ( warning_level(1:7) == 'WARNING' .or. &
              warning_level(1:7) == 'warning' ) then
            call mpp_set_warn_level ( WARNING )
    else
            call error_mesg ( 'utilities_init',  &
            'invalid entry for namelist variable warning_level', FATAL )
    endif

!--- open logfile (at beginning) and write version info ---

    unit = open_file ('logfile.out', action='write')
    if ( get_my_pe() == 0 ) then
         write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
         write (unit, nml=utilities_nml)
         write (unit,*) 'nml_error_codes=',  &
                         nml_error_codes(1:num_nml_error_codes)
    endif
    call close_file (unit)


end subroutine utilities_init

!#######################################################################

subroutine utilities_end

    call mpp_io_exit
    call mpp_domains_exit
    call mpp_exit
    call check_system_clock ('END OF RUN')

end subroutine utilities_end

!#######################################################################

subroutine close_file (unit, status)
   integer,          intent(in)           :: unit
   character(len=*), intent(in), optional :: status

   if (present(status)) then
       if (status(1:6) == 'delete' .or. status(1:6) == 'DELETE') then
           call mpp_close (unit, action=MPP_DELETE)
       else
           call mpp_close (unit)
       endif
   else
           call mpp_close (unit)
   endif
   

end subroutine close_file

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

end module utilities_mod

