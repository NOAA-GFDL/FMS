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

program test_affinity

!--- FMS modules
 use mpp_mod,          only: input_nml_file, mpp_error, mpp_pe, mpp_root_pe, FATAl
 use fms_affinity_mod

!--- namelist parameters
 character(len=15) :: component='COMPONENT      '
 integer:: nthreads   = 2
 logical:: use_hyper_thread = .false.
 namelist /test_affinity_nml/ component, nthreads, use_hyper_thread

!--- program vars
 integer:: io_status
 integer:: conc_threads
 character(len=32):: h_name

!-----------------------------------------------------------------------
!--- print initial test message
    print *, ''
    print *, '*** Testing affinity placement within FMS library...'
!--- initialize affinity
    call FMS_affinity_init()

    read(input_nml_file,test_affinity_nml, iostat=io_status)
    if (io_status > 0) then
      call mpp_error(FATAL, '=> test_affinity: Error reading input.nml')
    endif

    call fms_affinity_set (component, use_hyper_thread, nthreads)

!--- print success or failure message
    if (mpp_pe() == mpp_root_pe()) print *, '*** SUCCESS!'

end program test_affinity
