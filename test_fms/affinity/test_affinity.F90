!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

program test_affinity

!--- FMS modules
 use mpp_mod,          only: input_nml_file, mpp_error, mpp_pe, mpp_root_pe, FATAl
 use fms_mod,          only: fms_init, fms_end
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

     call fms_init()

!-----------------------------------------------------------------------
!--- print initial test message
    print *, ''
    print *, '*** Testing affinity placement within FMS library...'

    read(input_nml_file,test_affinity_nml, iostat=io_status)
    if (io_status > 0) then
      call mpp_error(FATAL, '=> test_affinity: Error reading input.nml')
    endif

    call fms_affinity_set (component, use_hyper_thread, nthreads)

!--- print success or failure message
    if (mpp_pe() == mpp_root_pe()) print *, '*** SUCCESS!'

    call fms_end()

end program test_affinity
