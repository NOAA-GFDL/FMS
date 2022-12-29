# Intel C
set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -sox -traceback")

set( CMAKE_C_FLAGS_RELEASE "-qno-opt-dynamic-align -O2 -debug minimal -fp-model source")

set( CMAKE_C_FLAGS_DEBUG "-O0 -g -ftrapuv")

set( CMAKE_C_LINK_FLAGS "")
