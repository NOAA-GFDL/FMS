program test_affinity

!--- FMS modules
 use mpp_mod,          only: input_nml_file, mpp_error, mpp_pe, mpp_root_pe, FATAl
 use fms_affinity_mod

!--- namelist parameters
 character(len=15) :: component='COMPONENT      '
 integer:: nthreads   = 2
 logical:: use_hyper_threads = .false.
 logical:: do_concurrent = .false.
 integer:: nthreads_2 = 2
 namelist /test_affinity_nml/ component, nthreads, use_hyper_threads, do_concurrent, nthreads_2

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

    conc_threads = nthreads
    if (do_concurrent) conc_threads = nthreads + nthreads_2

    call fms_set_affinity (component, conc_threads, use_hyper_threads)

!--- print success or failure message
    if (mpp_pe() == mpp_root_pe()) print *, '*** SUCCESS!'
 
end program test_affinity
