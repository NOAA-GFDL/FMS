# Intel LLVM based C

# flags from mkmf intel oneapi template
# Release is prod/opt flags, debug and repro are the same
set( CMAKE_C_FLAGS_BASE "-sox -traceback -march=core-avx-i")

set( isa_flags "-march=core-avx-i")

set( CMAKE_C_FLAGS "")

set( CMAKE_C_FLAGS_RELEASE "-O2 -debug minimal ${CMAKE_C_FLAGS_BASE} ${isa_flags}")
set( CMAKE_C_FLAGS_DEBUG "-O0 -g ${CMAKE_C_FLAGS_BASE}")
set( CMAKE_C_FLAGS_REPRO "-O2 -debug minimal ${CMAKE_C_FLAGS_BASE} ${isa_flags}")
set( CMAKE_C_FLAGS_NOFLAGS "")

set( CMAKE_C_LINK_FLAGS "")

# sets of flags used by ufs
set( CMAKE_C_FLAGS_DEBUGUFS "-O0 -g -ftrapv -traceback")
set( CMAKE_C_FLAGS_RELEASEUFS "-qno-opt-dynamic-align -O2 -debug minimal")

