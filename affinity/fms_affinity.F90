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
module fms_affinity_mod

  !--- standard system modules
  use iso_c_binding, only: c_int, c_bool
  use omp_lib

  !--- FMS modules
  use mpp_mod,    only: input_nml_file, mpp_pe, stdlog
  use fms_mod,    only: fms_init, check_nml_error, write_version_number, &
                        error_mesg, FATAL

  !--- default scoping
  implicit none
  private

  !--- namelist parameters
  logical:: debug_affinity = .false.
  logical(c_bool):: debug_cpuset = .false.
  namelist /fms_affinity_nml/ debug_affinity, debug_cpuset

  public fms_affinity_init, fms_set_affinity

  !---- version number
  ! Include variable "version" to be written to log file.
#include <file_version.h>

  logical :: module_is_initialized = .FALSE.

contains


  !--- initialization routine for affinity handling
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

    

  !--- routine to set affinity 
  subroutine fms_set_affinity (component, nthreads, use_hyper_threads)
    !--- interface variables
    character(len=*), intent(in):: component
    integer,          intent(in):: nthreads
    logical,          intent(in):: use_hyper_threads

    !--- local declarations for Fortran/C affinity interoperability
    integer(c_int):: get_cpuset_affinity
    integer(c_int):: get_cpu_affinity
    integer(c_int):: set_cpu_affinity
    integer(c_int):: cpuset_sz
    integer(c_int), dimension(:), allocatable:: cpu_set
    integer(c_int):: retcode

    !--- local variables
    integer:: th_num
    character(len=32):: h_name

     h_name = 'generic                         '

     if (use_hyper_threads) then
       cpuset_sz = nthreads
     else
       cpuset_sz = nthreads * 2
     endif
     allocate (cpu_set(0:cpuset_sz-1))
     retcode = get_cpuset_affinity(cpuset_sz, cpu_set, mpp_pe(), debug_cpuset)
     if (retcode == -1) then
       call error_mesg('fms_set_affinity',trim(component)//' fms_set_affinity: cpu_set size > allocated storage',FATAL)
     elseif (retcode < cpuset_sz) then
       call error_mesg('fms_set_affinity',trim(component)//' fms_set_affinity: cpu_set size smaller than expected',FATAL)
     endif
!$   call omp_set_num_threads(nthreads)
!$OMP PARALLEL NUM_THREADS(nthreads), PRIVATE(th_num, retcode)
!$   th_num = omp_get_thread_num() 
!$   retcode = set_cpu_affinity(cpu_set(th_num))
!$   if (retcode == -1) then
!$     call error_mesg('fms_set_affinigy',trim(component)//': issue setting cpu affinity', FATAL)
!$   endif
!$   if (debug_affinity) then
!$      call hostnm(h_name)
!$      print *, 'NOTE:',mpp_pe(),trim(component),' ',trim(h_name),get_cpu_affinity(),cpu_set(0),th_num
!$   endif
!$OMP END PARALLEL

  end subroutine fms_set_affinity

end module fms_affinity_mod
