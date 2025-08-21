# Precision-based Fortran compiler flags
set(r4_flags "-real-size 32") # Fortran flags for 32BIT precision
set(r8_flags "-real-size 64") # Fortran flags for 64BIT precision

# Minimal set of flags for stand release and debug build types
set(CMAKE_Fortran_FLAGS "${RELEASE} -safe-cray-ptr")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0")

# ufs flags to reproduce past behavior
set(ufs_flags_base "${CMAKE_Fortran_FLAGS} -fpp -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -align array64byte -nowarn -sox -traceback")

set(CMAKE_Fortran_FLAGS_RELEASE "${ufs_flags_base} -O2 -debug minimal -fp-model source -nowarn -qoverride-limits -qno-opt-dynamic-align -qopt-prefetch=3")
set(CMAKE_Fortran_FLAGS_DEBUGUFS "${ufs_flags_base} -g -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fpe0 -ftrapuv")

set(CMAKE_Fortran_LINK_FLAGS "")
