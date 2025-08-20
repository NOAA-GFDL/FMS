# Precision-based Fortran compiler flags
set(r4_flags "-real-size 32") # Fortran flags for 32BIT precision
set(r8_flags "-real-size 64") # Fortran flags for 64BIT precision

set(isa_flag "-march=core-avx-i")

# base flags from mkmf oneapi template
set(CMAKE_Fortran_FLAGS_BASE "-fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -i4 -nowarn -traceback")

# set the default as empty so we can add on to it
set(CMAKE_Fortran_FLAGS "")

set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -debug minimal -fp-model source ${isa_flag} ${CMAKE_Fortran_FLAGS_BASE}")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check -check noarg_temp_created -check nopointer -check nouninit -warn -warn noerrors -fpe0 -ftrapuv ${isa_flag} ${CMAKE_Fortran_FLAGS_BASE}")
set(CMAKE_Fortran_FLAGS_REPRO "-O2 -debug minimal -fp-model source ${isa_flag} ${CMAKE_Fortran_FLAGS_BASE}")
set(CMAKE_Fortran_FLAGS_NOFLAGS "")

set(CMAKE_Fortran_LINK_FLAGS "-fuse-ld=lld")

# ufs flags
set(CMAKE_Fortran_SHARED_FLAGS "-fpp -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -align array64byte -nowarn -traceback")
set(CMAKE_Fortran_FLAGS_RELEASEUFS "-O2 -debug minimal -nowarn -qoverride-limits -qno-opt-dynamic-align ${CMAKE_Fortran_SHARED_FLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUGUFS "-g -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fpe0 -ftrapuv ${CMAKE_Fortran_SHARED_FLAGS}")
