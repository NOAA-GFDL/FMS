!-----------------------------------------------------------------------
!                  CLOCKS module (for timing code sections)
!
! AUTHOR: V. Balaji (vb@gfdl.gov)
!         SGI/GFDL Princeton University
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.  
!-----------------------------------------------------------------------

!removed system_clock_sgi overload: the hires clock on SGI (800ns) overflows 
!       a 4-byte integer tick very quickly (2^32 ticks is only about 3200 s).
!mpp_clocks has alternate formulation, see mpp.F90
!#ifdef __sgi
!#define SYSTEM_CLOCK system_clock_sgi
!#endif

module clocks_mod
! <CONTACT EMAIL="vb@gfdl.noaa.gov">
!       V. Balaji
! </CONTACT>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!    <TT>clocks_mod</TT> is a set of simple calls for timing f90
!    code and code sections.
! </OVERVIEW>
! <DESCRIPTION>
!    In parallel environments, the key timing information is the
!    wallclock ("real") time, not the CPU time, of a run. F90 provides the
!    <TT>system_clock(3F)</TT> intrinsic to retrieve timing
!    information from the system realtime clock.
!
!    <TT>clocks_mod</TT> uses the F90
!    <TT>system_clock(3F)</TT> intrinsic to measure time. The main
!    call is to the function <TT>tick()</TT>. Clocks can be set up
!    to provide direct timing of a section of code, or for cumulative
!    timing of code sections within loops.
!
!    The overhead of calls to the system clock is typically measured in
!    microseconds. However, the resolution of the clock may be higher or
!    lower than this overhead. The resolution is printed when the module is
!    initialized. A test program is supplied with the module which, among
!    other things, measures the calling overhead.
!
!    On SGI systems <TT>SYSTEM_CLOCK</TT> is transparently
!    overloaded with a higher resolution clock made available in a
!    non-portable fortran interface made available by
!    <TT>nsclock.c</TT>. This approach will eventually be extended to other
!    platforms.
!
!    This module has now been extended to work in parallel environments,
!    using the <LINK SRC="http://www.gfdl.gov/~vb/mpp.html"><TT>mpp</TT>
!    package</LINK>.  In a parallel environment, the clocks are synchronized
!    across all the PEs in the current pelist at the top of the timed code
!    section, but allows each PE to complete the code section at different
!    times. This allows us to measure load imbalance for a given code
!    section. Statistics are written to <TT>stdout</TT> by
!    <TT>clocks_exit</TT>.
!
!    While the nesting of clocks is allowed, please note that
!    synchronization on an inner clock may distort outer clock measurements
!    of load imbalance.
! </DESCRIPTION>

  use mpp_mod, only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe, mpp_error, mpp_sync
  use mpp_mod, only: mpp_max, mpp_min, mpp_sum
  use mpp_mod, only: stdlog, stdout, stderr
  use mpp_mod, only: WARNING, FATAL

  use platform_mod, only: i8_kind

  use fms_mod, only: write_version_number

  implicit none
  private
  integer, private :: ticks_per_sec, max_ticks, ref_tick, start_tick, end_tick
  real, private :: tick_rate
  integer, private, parameter :: MAX_CLOCKS=256
  integer, private :: clock_num=0
  logical, private :: module_is_initialized=.FALSE., verbose=.FALSE.
  integer, parameter, public :: CLOCKS_VERBOSE=1
  character(len=128), private :: errortxt
!clocks are stored in this internal type
  type, private :: clock
     character(len=24) :: name
     integer :: ticks, calls
  end type clock
  type(clock) :: clocks(MAX_CLOCKS)
  integer :: id_cumul_clock

  public :: clocks_init, clocks_exit, get_clock, clock_id, tick

  character(len=128), private :: &
   version='$Id: clocks.F90,v 2.5 2003/04/09 21:15:38 fms Exp $'
  character(len=128), private :: &
   tagname='$Name: inchon $'

  contains

#ifdef __sgi
    subroutine system_clock_sgi( count, count_rate, count_max )
!mimics F90 SYSTEM_CLOCK intrinsic
      integer, intent(out), optional :: count, count_rate, count_max
      integer(kind=i8_kind) :: sgi_tick, sgi_ticks_per_sec, sgi_max_tick
      if( PRESENT(count) )then
          count = sgi_tick()
      end if
      if( PRESENT(count_rate) )then
          count_rate = sgi_ticks_per_sec()
      end if
      if( PRESENT(count_max) )then
          count_max = 2**sgi_max_tick()
      end if
      return
    end subroutine system_clock_sgi
#endif

! <SUBROUTINE NAME="clocks_init">

!   <OVERVIEW>
!     Initialize clocks module.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Called to initialize the <TT>clocks_mod</TT> package. Some
!     information is printed regarding the version of this module and the
!     resolution of the system clock. <TT>flag</TT> may be used, for
!     example, in a parallel run, to have only one of the PEs print this
!     information.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call clocks_init(flag)
!   </TEMPLATE>

!   <IN NAME="flag" TYPE="integer">
!     if flag is set, only print if flag=0
!     for instance, flag could be set to pe number by the calling program
!     to have only PE 0 in a parallel run print info
!   </IN>

    subroutine clocks_init(flags)
!flag is used to set verbose flag, currently unused
      integer, intent(in), optional :: flags
      integer :: i

      if( module_is_initialized )return
      module_is_initialized = .TRUE.

      if( PRESENT(flags) )then
          verbose = flags.EQ.CLOCKS_VERBOSE
      end if
      call mpp_init()
      call mpp_sync()
!initialize clocks and reference tick
      call SYSTEM_CLOCK( start_tick, ticks_per_sec, max_ticks )
      tick_rate = 1./ticks_per_sec
      call write_version_number( version, tagname )
      if( mpp_pe().EQ.mpp_root_pe() )then
!write version
         write( stdlog(), '(a,es12.4,a,i10,a)' ) &
              'Realtime clock resolution=', tick_rate, ' sec (', ticks_per_sec, ' ticks/sec)'
      end if
!default clock name is Clock001, etc
      do i = 1,MAX_CLOCKS
         write( clocks(i)%name,'(a5,i3.3)' ) 'Clock', i
         clocks(i)%ticks = 0
         clocks(i)%calls = 0
      end do
      id_cumul_clock = clock_id( 'Total measured time' )
      ref_tick = start_tick

      return
    end subroutine clocks_init
! </SUBROUTINE>

! <FUNCTION NAME="clock_id">

!   <OVERVIEW>
!     Return an ID to a new or existing cumulative clock. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     This is used to return an ID to a clock that may be used for timing a
!     code section. Currently up to 256 (an arbitrarily chosen setting for the
!     internal parameter <TT>max_clocks</TT>) clocks can be set. The cumulative
!     times can be printed at the end of the run by a call to
!     <LINK SRC="#clocks_exit"><TT>clocks_exit().</TT></LINK> The name can be a 
!     new or existing clock.
!     <I>Note that <TT>name</TT> is restricted to 24 characters</I>. If you 
!     enter a longer
!     name, it is silently and gracefully truncated. This can be problematic if
!     you inadvertently give different clocks names that differ only beyond the
!     24th character. These will look to the clocks module as the same clock.
!   </DESCRIPTION>
!   <TEMPLATE> clock_id(name) </TEMPLATE>

!   <IN NAME="name" TYPE="character(len=*)"> </IN>

    function clock_id(name)
      integer :: clock_id
      character(len=*), intent(in) :: name

      if( .NOT.module_is_initialized ) &
           call mpp_error( FATAL, 'CLOCKS: must first call clocks_init().' )
      clock_id = 1
      do while( trim(name).NE.trim(clocks(clock_id)%name) )
         if( clock_id.GT.clock_num )then
             if( clock_num.LT.MAX_CLOCKS )then
                 clock_num = clock_id
                 clocks(clock_id)%name = name
             else
                 write( errortxt, '(a,i4)' )'CLOCKS: you are requesting too many clocks, max clocks=', MAX_CLOCKS
                 call mpp_error( WARNING, trim(errortxt) )
             end if
             return
         end if
         clock_id = clock_id + 1
      end do
    end function clock_id
! </FUNCTION>

! <FUNCTION NAME="tick">

!   <OVERVIEW>
!     Return time on system clock. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     <TT>tick</TT> returns the current tick of the system clock.
!     If <TT>string</TT> is present, it prints the time elapsed since the 
!     reference tick.<BR/>
!     Otherwise if <TT>id</TT> is present, the time elapsed since the reference tick is
!     accumulated to the clock <TT>id</TT>.<BR/>
!     Otherwise if <TT>name</TT> is present, the time elapsed since the reference tick
!     is accumulated to the clock whose name is <TT>name</TT>. There is slightly
!     larger overhead for this option, to resolve the name to an ID. It is
!     recommended to use <TT>id=</TT>, especially for small sections.<BR/>
!     The reference tick is either the value of the system clock at the last call
!     to <TT>tick()</TT> (or <TT>clocks_init()</TT>), or else as given by the optional
!     <TT>since</TT> argument.
!   </DESCRIPTION>
!   <TEMPLATE>
!     tick( string, id, name, since )
!   </TEMPLATE>

!   <IN NAME="string" TYPE="character(len=*)"></IN>
!   <IN NAME="name" TYPE="character(len=*)"></IN>
!   <IN NAME="id" TYPE="integer"></IN>
!   <IN NAME="since" TYPE="integer"></IN>

    function tick( string, id, name, since )
      integer :: tick
      character(len=*), intent(in), optional :: string
      character(len=*), intent(in), optional :: name
      integer, intent(in), optional :: id, since
      integer :: current_tick, cid

      if( .NOT.module_is_initialized ) &
           call mpp_error( FATAL, 'CLOCKS: must first call clocks_init().' )
!take time first, so that this routine's overhead isn't included
      call SYSTEM_CLOCK(current_tick)
!ref_tick is the clock value at the last call to tick (or clocks_init)
!unless superseded by the since argument.
      if( PRESENT(since) )ref_tick = since
!correct ref_tick in the unlikely event of clock rollover
      if( current_tick.LT.ref_tick )ref_tick = ref_tick - max_ticks
      if( PRESENT(string) )then
!print time since reference tick
          write( stdout(),'(a,f14.6)' ) &
               'CLOCKS: '//trim(string), (current_tick-ref_tick)*tick_rate
      else if( PRESENT(id) )then
          cid = id
!accumulate time on clock id
          if( 0.LT.cid .AND. cid.LE.MAX_CLOCKS )then
              clocks(cid)%ticks = clocks(cid)%ticks + current_tick - ref_tick
              clocks(cid)%calls = clocks(cid)%calls + 1
          else
              call mpp_error( FATAL, 'CLOCKS: invalid id' )
          end if
      else if( PRESENT(name) )then
          cid = clock_id(name)
!accumulate time on clock id
          if( 0.LT.cid .AND. cid.LE.MAX_CLOCKS )then
              clocks(cid)%ticks = clocks(cid)%ticks + current_tick - ref_tick
              clocks(cid)%calls = clocks(cid)%calls + 1
          else
              call mpp_error( FATAL, 'CLOCKS: invalid id' )
          end if
      end if
!reset reference tick
      call mpp_sync()
      call SYSTEM_CLOCK(ref_tick)
      tick = ref_tick

      return
    end function tick
! </FUNCTION>

! <SUBROUTINE NAME="get_clock">

!   <OVERVIEW>
!     Retrieve information from a cumulative clock.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This is used to return information stored on the clock whose ID is
!     <TT>id</TT>. The subroutine returns any or all of the
!     information held in the following: <TT>ticks</TT>, for the
!     total accumulated clock ticks for this ID; <TT>calls</TT>,
!     the number of intervals measured with this clock (i.e, the number of
!     times <TT>tick()</TT> was called with this ID);
!     <TT>total_time</TT>, the total time in seconds on this clock;
!     and <TT>time_per_call</TT>, the time per measured interval on
!     this clock (<TT>total_time/calls</TT>). This routine is used if you wish to retrieve this information in a
!     variable. Otherwise, <TT>clocks_exit()</TT> may be used to
!     print this information at termination.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_clock( id, ticks, calls, total_time, time_per_call)
!   </TEMPLATE>

!   <IN NAME="id" TYPE="integer"> </IN>
!   <OUT NAME="ticks" TYPE="integer"></OUT>
!   <OUT NAME="calls" TYPE="integer"></OUT>
!   <OUT NAME="total_time" TYPE="real"></OUT>
!   <OUT NAME="time_per_call" TYPE="real"></OUT>

    subroutine get_clock( id, ticks, calls, total_time, time_per_call )
      integer, intent(in) :: id
      integer, intent(out), optional :: ticks, calls
      real, intent(out), optional :: total_time, time_per_call

      if( .NOT.module_is_initialized ) &
           call mpp_error( FATAL, 'CLOCKS: must first call clocks_init().' )
      if( 0.LT.id .AND. id.LE.MAX_CLOCKS )then
          if( PRESENT(ticks) )ticks = clocks(id)%ticks
          if( PRESENT(calls) )calls = clocks(id)%calls
          if( PRESENT(total_time) )total_time = clocks(id)%ticks*tick_rate
          if( PRESENT(time_per_call) )time_per_call = &
               clocks(id)%ticks*tick_rate/clocks(id)%calls
      else
          call mpp_error( FATAL, 'CLOCKS: invalid id' )
      end if

      return
    end subroutine get_clock
! </SUBROUTINE>

    subroutine clocks_exit()

! <SUBROUTINE NAME="clocks_exit">

!   <OVERVIEW>
!     Exit clocks_mod.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This prints the values of all cumulative clocks. In a parallel
!     environment, statistics across PEs are also printed.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call clocks_exit()
!   </TEMPLATE>
      integer :: i
      real :: t, tmax, tmin, tavg, tstd, tavg_call
      real :: t_total

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'CLOCKS: must first call clocks_init.' )

      call mpp_sync()
      end_tick = tick( id=id_cumul_clock, since=start_tick )
      if( mpp_pe().EQ.mpp_root_pe() )write( stdout(),'(32x,a)' ) &
           '   calls          tmin          tmax          tavg          tstd     tavg/call  tfrac'
      do i = 1, clock_num
         t = clocks(i)%ticks*tick_rate
         if( clocks(i)%ticks.LT.0 )then
             write( stderr(),'(a,i4,2i16)' )'pe, clocks(i)%ticks=', mpp_pe(), clocks(i)%ticks
             call mpp_error( FATAL, 'CLOCKS_EXIT: negative total_time.' )
         end if
         tmax = t; call mpp_max(tmax)
         tmin = t; call mpp_min(tmin)
         tavg = t;           call mpp_sum(tavg); tavg = tavg/mpp_npes()
         tstd = (t-tavg)**2; call mpp_sum(tstd); tstd = sqrt(tstd/mpp_npes())
         if( i.EQ.1 )t_total = tavg
         if( mpp_pe().EQ.mpp_root_pe() )then
             write( stdout(),'(a32,i8,5f14.6,f7.3)' ) &
                  'CLOCKS: '//clocks(i)%name, &
                  clocks(i)%calls, tmin, tmax, tavg, tstd, tavg/max(clocks(i)%calls,1), tavg/t_total
         end if
      end do

      return
    end subroutine clocks_exit
! </SUBROUTINE>

end module clocks_mod

#ifdef test_clocks
#ifdef SYSTEM_CLOCK
#undef SYSTEM_CLOCK
#endif
program test
  use clocks_mod
  use mpp_mod, only: mpp_init, mpp_exit, stdout

  call mpp_init( log=stdout() )
  call clocks_init()
  k = tick()
  call pxfsleep(1,i,j)
  k = tick( 'Time 1 second', since=k )

  call SYSTEM_CLOCK(count_rate=irate)
  call SYSTEM_CLOCK(k)
  call pxfsleep(1,i,j)
  call SYSTEM_CLOCK(j)
  print *, 'Direct time ', float(j-k)/irate

  id = clock_id( 'Clock1' )
  k = tick()
  call pxfsleep(1,i,j)
  j = tick( id=id, since=k )
  call get_clock( id, total_time=t )
  print *, 't=', t
  
  k = tick()
  do n = 1,1000000
     j = tick()
  end do

  j = tick( 'Time per 1000000 calls', since=k )
  call clocks_exit()
  call mpp_exit
end program test
#endif

! <INFO>

!   <COMPILER NAME="COMPILING AND LINKING SOURCE">
!    Any module or program unit using <TT>clocks_mod</TT> must contain the line
!
!    <PRE>use clocks_mod</PRE>
!
!
!   The parallel version of this module requires the <LINK
!   SRC="http://www.gfdl.gov/~vb/mpp.html"><TT>mpp</TT> package</LINK>.
!
!   Compiling with the cpp flag <TT>test_clocks</TT> turned on:
!
!    <PRE>f90 -Dtest_clocks clocks.F90</PRE>
!
!   will produce a program that will exercise certain portions of the
!   <TT>clocks_mod</TT> module.
!   </COMPILER>
!   <PRECOMP FLAG="PORTABILITY">      
!     <TT>clocks_mod</TT> is fully f90 standard-compliant. There are
!     no portability issues.
!
!     On SGI systems, the <TT>f90</TT> standard <TT>SYSTEM_CLOCK</TT>
!     intrinsic is overloaded with a non-portable fortran interface to a
!     higher-precision clock. This is distributed with the MPP package as
!     <TT>nsclock.c</TT>. This approach will eventually be extended to other
!     platforms, since the resolution of the default clock is often too
!     coarse for our needs.
!   </PRECOMP>
!   <LOADER FLAG="">
!       ACQUIRING SOURCE<\BR>
!
!       GFDL users can check it out of the main CVS repository as the
!       <TT>clocks</TT> CVS module. The current public tag is <TT>fez</TT>.
!       External users can download the source <LINK
!       SRC="ftp://ftp.gfdl.noaa.gov/pub/vb/utils/clocks.F90">here</LINK>.
!   </LOADER>
!   <TESTPROGRAM NAME="">  test </TESTPROGRAM>
!   <BUG>
!     The <TT>SYSTEM_CLOCK</TT> intrinsic has a limited range before the
!     clock rolls over. The maximum time interval that may be measured
!     before rollover depends on the default integer precision, and is
!     <TT>COUNT_MAX/COUNT_RATE</TT> seconds. Timing a code section longer
!     than this interval will give incorrect results. The <TT>clocks</TT>
!     entry in the logfile reports the rollover time interval. Note that
!     this is a limitation, or "feature" of the <TT>f90 SYSTEM_CLOCK</TT>
!     intrinsic.
!   </BUG>
!   <NOTE>
!
!     <LINK SRC="clocks.F90">Using <TT>clocks_mod</TT></LINK>
!
!   
!     In the simplest method of calling, just designate sections of a main
!     program with calls to <TT>tick()</TT>. <TT>tick()</TT>
!     by default measures time since the last call to
!     <TT>tick()</TT> (or to <TT>clocks_init()</TT>).
!
!     <PRE>
!     program main
!     call clocks_init()
!     !code section 1
!     ...
!     i = tick( 'code section 1' )
!     !code section 2
!     ...
!     i = tick( 'code section 2' )
!     !code section 3
!     ...
!     i = tick( 'code section 3' )
!     end
!     </PRE>
!
!     This will return timing information for the three regions of the
!     main program.
!
!     If, however, a subroutine in one of the code sections itself
!     contained calls to <TT>tick()</TT>, this would produce
!     erroneous information (since "the last call to <TT>tick()</TT>"
!     might refer to a call elsewhere). In this case, we set the reference
!     tick using the <TT>since</TT> argument to
!     <TT>tick()</TT>:
!
!     <PRE>
!     program main
!     call clocks_init()
!     i = tick()
!     !code section 1
!     ...
!     i = tick( 'code section 1', since=i )
!     !code section 2
!     ...
!     i = tick( 'code section 2', since=i )
!     !code section 3
!     ...
!     i = tick( 'code section 3', since=i )
!     end
!     </PRE>
!
!     A third way to use <TT>clocks_mod</TT> is to produce cumulative
!     times of code sections within loops. Here we first call
!     <TT>clock_id</TT> to set up a clock with an
!     <TT>id</TT>, and accumulate times to this ID.
!
!     <PRE>
!     program main
!     call clocks_init()
!     id1 = clock_id( 'Code section 1' )
!     id2 = clock_id( 'Code section 2' )
!     id3 = clock_id( 'Code section 3' )
!     do j = 1,10000
!        i = tick()
!     !code section 1
!     ...
!        i = tick( id=id1, since=i )
!     !code section 2
!     ...
!        i = tick( id=id2, since=i )
!     !code section 3
!     ...
!        i = tick( id=id3, since=i )
!     end do
!     call clocks_exit()
!     end
!     </PRE>
!
!     The call to <TT>clocks_exit</TT> above prints the timings
!     for the code sections.

!   </NOTE>

! </INFO>

