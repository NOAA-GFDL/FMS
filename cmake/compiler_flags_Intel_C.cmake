# Intel C
set( CMAKE_C_FLAGS "")

set( base_flags "-sox -traceback")

set( isa_flags "-march=core-avx-i -qno-opt-dynamic-align")

set( CMAKE_C_FLAGS_RELEASE "-O2 -debug minimal ${base_flags} ${isa_flags}")
set( CMAKE_C_FLAGS_REPRO "-O2 -debug minimal ${base_flags} ${isa_flags}")
set( CMAKE_C_FLAGS_DEBUG "-O0 -g -ftrapuv ${base_flags} ${isa_flags}")

set( CMAKE_C_LINK_FLAGS "")

# ufs flags
set( CMAKE_C_FLAGS_RELEASEUFS "-qno-opt-dynamic-align -O2 -debug minimal -fp-model source")
set( CMAKE_C_FLAGS_DEBUGUFS "-O0 -g -ftrapuv")
