# Precision-based Fortran compiler flags
set(r4_flags "-real-size 32") # Fortran flags for 32BIT precision
set(r8_flags "-real-size 64") # Fortran flags for 64BIT precision
set(r8_flags "${r8_flags} -no-prec-div -no-prec-sqrt")

# Intel Fortan
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -align array64byte -nowarn -sox -traceback")

set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -debug minimal -fp-model source -nowarn -qoverride-limits -qno-opt-dynamic-align -qopt-prefetch=3")

set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fpe0 -ftrapuv")

set(CMAKE_Fortran_LINK_FLAGS "")
