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

!> @defgroup fms_affinity_mod fms_affinity_mod
!> @ingroup affinity
!> @brief Fortran API interfaces to set the thread affinity.
!! API interfaces to allow setting and getting thread affinity.  The thread affinity get and set
!! are managed in the C routines in affinity.c.
!!
!! @author Rusty Benson

!> @file
!> File for @ref fms_affinity_mod

!> @addtogroup fms_affinity_mod
!> @{
module fms_affinity_mod
  !--- standard system modules
  use, intrinsic :: iso_c_binding, only: c_int, c_bool
  use omp_lib

  !--- FMS modules
  use mpp_mod,    only: input_nml_file, mpp_pe, stdlog
  use fms_mod,    only: fms_init, check_nml_error, write_version_number, &
                        error_mesg, FATAL, NOTE

  !--- default scoping
  implicit none
  private

  interface

    !> Interface to get affinity from the current component.
    !!
    !> Defined in @ref affinity.c.
    !> @ingroup fms_affinity_mod
    integer(KIND=c_int) function fms_affinity_get() bind(c, name="get_cpu_affinity")
      import c_int
    end function fms_affinity_get

    !> Private interface to retrieve this groups CPU set and it's size.
    !!
    !> Defined in @ref affinity.c.
    integer(KIND=c_int) function get_cpuset(fsz, output, pe, debug) bind(c, name="get_cpuset")
      import c_int, c_bool
      integer(KIND=c_int), value, intent(in) :: fsz, pe
      integer(KIND=c_int), dimension(*), intent(inout) :: output
      logical(KIND=c_bool), value :: debug
    end function get_cpuset

    !> Private interface to set CPU afinity to a given core.
    !!
    !> Defined in @ref affinity.c.
    integer(KIND=c_int) function set_cpu_affinity(cpu) bind(c, name="set_cpu_affinity")
      import c_int
      integer(KIND=c_int), value, intent(in) :: cpu
    end function set_cpu_affinity
  end interface

  !--- namelist parameters
  logical:: affinity = .true.
  logical:: strict = .true.
  logical:: debug_affinity = .false.
  logical(c_bool):: debug_cpuset = .false.
  namelist /fms_affinity_nml/ affinity, strict, debug_affinity, debug_cpuset

  public fms_affinity_init, fms_affinity_get, fms_affinity_set

  !---- version number
  ! Include variable "version" to be written to log file.
#include <file_version.h>

  logical :: module_is_initialized = .FALSE.

contains

  !> Initialization routine for affinity handling
  subroutine fms_affinity_init()
    !--- local variables
    integer:: io_stat
    integer:: ierr
    integer:: unit

    !--- return if module is initialized
    if (module_is_initialized) return

    !--- ensure fms/mpp has been initialized
    call fms_init()

    !--- read in namelist
    read(input_nml_file, fms_affinity_nml, iostat=io_stat)
    ierr = check_nml_error(io_stat,'fms_affinity_nml')

    !--- output information to logfile
    call write_version_number("fms_affinity_mod", version)
    unit = stdlog()
    write(unit,nml=fms_affinity_nml)

    module_is_initialized = .TRUE.

  end subroutine fms_affinity_init


  !> Routine to set affinity for a component
  subroutine fms_affinity_set (component, use_hyper_thread, nthreads)
    !--- interface variables
    character(len=*),  intent(in):: component !< Component name
    logical,           intent(in):: use_hyper_thread !< .TRUE. if using hyperthreads
    integer,           intent(in):: nthreads !< Number of threads

    !--- local declarations for Fortran/C affinity interoperability
    integer(c_int):: cpuset_sz
    integer(c_int), dimension(:), allocatable:: cpu_set
    integer(c_int):: retcode

    !--- local variables
    character(len=32):: h_name
    integer:: MSG_TYPE
    integer:: th_num
    integer:: indx

    if (.not. affinity) return

    if (strict) then
      MSG_TYPE = FATAL
    else
      MSG_TYPE = NOTE
    endif

    h_name = 'generic                         '

    !--- allocate storage for cpuset
    if (use_hyper_thread) then
      cpuset_sz = nthreads
    else
      cpuset_sz = nthreads * 2
    endif
    allocate (cpu_set(0:cpuset_sz-1))

    !--- get cpuset for this MPI-rank
    retcode = get_cpuset(cpuset_sz, cpu_set, mpp_pe(), debug_cpuset)
    if (retcode == -1) then
      call error_mesg('fms_affinity_set',trim(component)//' cpu_set size > allocated storage',FATAL)
    elseif ( (retcode == cpuset_sz/2) .and. (retcode == nthreads) ) then
      call error_mesg('fms_affinity_set',trim(component)//' affinity assumes hyper-threading hardware disabled',NOTE)
    elseif (retcode < cpuset_sz) then
      call error_mesg('fms_affinity_set',trim(component)//' cpu_set size smaller than expected',MSG_TYPE)
    endif

    !--- set the affinity for the MPI-rank
    retcode = set_cpu_affinity(cpu_set(0))
    if (retcode == -1) then
      call error_mesg('fms_affinity_set',trim(component)//': issue setting cpu affinity', FATAL)
    endif

    !--- set affinity for threads associated with this MPI-rank
!$OMP PARALLEL NUM_THREADS (nthreads) &
!$OMP&         DEFAULT (none) &
!$OMP&         SHARED (use_hyper_thread, cpuset_sz, component, cpu_set, debug_affinity) &
!$OMP&         PRIVATE (th_num, indx, retcode, h_name)
!$  th_num = omp_get_thread_num()
    !--- handle hyper threading case by alternating threads between logical and virtual cores
!$  if (use_hyper_thread) then
!$    if (mod(th_num,2) == 0 ) then
!$      indx = th_num/2
!$    else
!$      indx = (cpuset_sz - 1) - ((cpuset_sz - 1) - th_num)/2
!$    endif
!$  else
!$    indx = th_num
!$  endif
!$  retcode = set_cpu_affinity(cpu_set(indx))
!$  if (retcode == -1) then
!$    call error_mesg('fms_affinity_set',trim(component)//': issue setting cpu affinity', FATAL)
!$  endif
    !--- output affinity placement
!$  if (debug_affinity) then
!$    call hostnm(h_name)
!$    print *, 'DEBUG:',mpp_pe(),trim(component),' ',trim(h_name),fms_affinity_get(),th_num
!$  endif
!$OMP END PARALLEL

  end subroutine fms_affinity_set
end module fms_affinity_mod
