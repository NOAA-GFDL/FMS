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
   version='$Id: clocks.F90,v 2.4 2002/07/16 22:54:36 fms Exp $'
  character(len=128), private :: &
   tagname='$Name: havana $'

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

    subroutine clocks_init(flags)
!initialize clocks module
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

    function clock_id(name)
!return an ID for a new or existing clock
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

    subroutine clocks_exit()
!print all clocks
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
