project(fft Fortran)

set(${PROJECT_NAME}_SOURCES
  fft.F90
  )
set(${PROJECT_NAME}99_SOURCES
  fft99.F90
  )

# The target ${PROJECT_NAME}99 is needed due to an error in how cmake
# calculates dependencies in Fortran.  It cannot correctly get the
# dependency on something like:
# #ifndef SGICRAY
# #ifndef NAGFFT
# use    fft99_mod, only: fft991, set99
# #endif
# #endif
add_library(${PROJECT_NAME}99 OBJECT ${${PROJECT_NAME}99_SOURCES})
add_dependencies(${PROJECT_NAME}99 constants mpp)
add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants platform mpp fms_mod fft99)

