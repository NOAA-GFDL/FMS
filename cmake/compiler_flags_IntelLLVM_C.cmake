# Intel LLVM based C

set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

set( CMAKE_C_FLAGS_RELEASE "-O2")
set( CMAKE_C_FLAGS_DEBUG "-O0 -g -traceback")

set( CMAKE_C_LINK_FLAGS "")

# sets of flags used by ufs
set( CMAKE_C_FLAGS_DEBUGUFS "-O0 -g -ftrapv -traceback")
set( CMAKE_C_FLAGS_RELEASEUFS "-qno-opt-dynamic-align -O2 -debug minimal")

