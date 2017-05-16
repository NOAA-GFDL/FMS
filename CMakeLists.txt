# -*-cmake-*-

# Minimul version of CMake required
cmake_minimum_required (VERSION 3.7)

# Indicate where the Fortran module files should be placed during the
# build.
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

project(libfms C Fortran)

include(${CMAKE_ROOT}/Modules/CheckFortranCompilerFlag.cmake)
include(${CMAKE_ROOT}/Modules/CheckCCompilerFlag.cmake)

# Default values for variables
set(CMAKE_MODULE_PATH ${libfms_SOURCE_DIR}/CMakeModules ${CMAKE_MODULE_PATH})
set(OPENMP TRUE CACHE BOOL "To determine if OpenMP build should be performed")
set(BUILD_TARGET "PROD" CACHE STRING "Build Target Type.  Can be \"PROD|REPRO|DEBUG\".  Use -DCMAKE_BUILD_TYPE=DEBUG for debug builds")

# Force CMAKE_BUILD_TYPE to DEBUG if BUILD_TARGET STREQUAL DEBUG
if (BUILD_TARGET STREQUAL DEBUG)
  set(CMAKE_BUILD_TYPE "DEBUG")
else (BUILD_TARGET STREQUAL DEBUG)
  set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "" FORCE)
endif (BUILD_TARGET STREQUAL DEBUG)

# Fortran compiler flags
if (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  set(CMAKE_Fortran_FLAGS
    "-fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -i4 -r8 -nowarn -g -sox -traceback")
  set(CMAKE_Fortran_FLAGS_PROD
    "-O2")
  set(CMAKE_Fortran_FLAGS_REPRO
    "-O2 -fltconsistency")
  set(CMAKE_Fortran_FLAGS_DEBUG
    "-O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -debug variable_locations -fpe0 -ftrapuv")
  set(CMAKE_Fortran_FLAGS_VERBOSE
    "-v -V -what -warn all")
  # Determine the OpenMP flag for the compiler
  # Newer versions of Intel use -qopenmp
  CHECK_Fortran_COMPILER_FLAG(-qopenmp _has_qopenmp)
  if (_has_qopenmp)
    set(CMAKE_Fortran_FLAGS_OPENMP
      "-qopenmp")
  else(_has_qopenmp)
    set(CMAKE_Fortran_FLAGS_OPENMP
      "-openmp")
  endif (_has_qopenmp)
  # CPP/FPP Macros
  set(CMAKE_Fortran_CPPFLAGS
    "-DINTERNAL_FILE_NML")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set(CMAKE_Fortran_FLAGS
    "-fcray-pointer -fdefault-double-8 -fdefault-real-8 -Waliasing -ffree-line-length-none -fno-range-check")
  set(CMAKE_Fortran_FLAGS_PROD
    "-O3")
  set(CMAKE_Fortran_FLAGS_REPRO
    "-O2 -fbounds-check")
  set(CMAKE_Fortran_FLAGS_DEBUG
    "-O0 -g -W -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow")
  set(CMAKE_Fortran_FLAGS_VERBOSE
    "-Wall -Wextra")
  set(CMAKE_Fortran_FLAGS_OPENMP
    "-fopenmp")
  # CPP/FPP Macros
  set(CMAKE_Fortran_CPPFLAGS
    "-DINTERNAL_FILE_NML")
endif (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

# C compiler flags
if (CMAKE_C_COMPILER_ID STREQUAL Intel)
  set(CMAKE_C_FLAGS
    "-D__IFC -sox -traceback")
  set(CMAKE_C_FLAGS_PROD
    "-O2")
  set(CMAKE_C_FLAGS_REPRO
    "")
  set(CMAKE_C_FLAGS_DEBUG
    "-O0 -g -ftrapuv")
  set(CMAKE_C_FLAGS_VERBOSE
    "-w3")
  # Determine the OpenMP flag for the compiler
  # Newer versions of Intel use -qopenmp
  CHECK_C_COMPILER_FLAG(-qopenmp _has_qopenmp)
  if (_has_qopenmp)
    set(CMAKE_C_FLAGS_OPENMP
      "-qopenmp")
  else(_has_qopenmp)
    set(CMAKE_C_FLAGS_OPENMP
      "-openmp")
  endif (_has_qopenmp)
elseif (CMAKE_C_COMPILER_ID STREQUAL GNU)
  set(CMAKE_C_FLAGS
    "-D__IFC")
  set(CMAKE_C_FLAGS_PROD
    "-O2")
  set(CMAKE_C_FLAGS_REPRO
    "")
  set(CMAKE_C_FLAGS_DEBUG
    "-O0 -g")
  set(CMAKE_C_FLAGS_VERBOSE
    "-Wall -Wextra")
  set(CMAKE_C_FLAGS_OPENMP
    "-fopenmp")
endif (CMAKE_C_COMPILER_ID STREQUAL Intel)

# OpenMP is on by default.
if (OPENMP)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_OPENMP}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_OPENMP}")
endif (OPENMP)

if (CMAKE_BUILD_TYPE STREQUAL RELEASE)
  set(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_${BUILD_TARGET}} CACHE STRING "" FORCE)
  set(CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_${BUILD_TARGET}} CACHE STRING "" FORCE)
endif (CMAKE_BUILD_TYPE STREQUAL RELEASE)

# Required external packages/libraries
find_package(NetCDF)
find_package(MPI)

# Include directory required by all
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${PROJECT_SOURCE_DIR}/include
  )

# **********************************************************************
# Source Files

add_subdirectory(constants)
add_subdirectory(platform)
add_subdirectory(tridiagonal)
add_subdirectory(mpp)
add_subdirectory(memutils)
add_subdirectory(drifters)
add_subdirectory(fms)
add_subdirectory(axis_utils)
add_subdirectory(block_control)
add_subdirectory(fft)
add_subdirectory(field_manager)
add_subdirectory(mosaic)
add_subdirectory(sat_vapor_pres)
add_subdirectory(time_manager)
add_subdirectory(tracer_manager)
add_subdirectory(astronomy)
add_subdirectory(column_diagnostics)
add_subdirectory(diag_manager)
add_subdirectory(horiz_interp)
add_subdirectory(oda_tools)
add_subdirectory(random_numbers)
add_subdirectory(time_interp)
add_subdirectory(topography)
add_subdirectory(amip_interp)
add_subdirectory(coupler)
add_subdirectory(data_override)
add_subdirectory(exchange)
add_subdirectory(interpolator)

add_library(fms
   $<TARGET_OBJECTS:constants>
   $<TARGET_OBJECTS:platform>
   $<TARGET_OBJECTS:tridiagonal>
   $<TARGET_OBJECTS:mpp>
   $<TARGET_OBJECTS:memutils>
   $<TARGET_OBJECTS:drifters>
   $<TARGET_OBJECTS:fms_mod>
   $<TARGET_OBJECTS:axis_utils>
   $<TARGET_OBJECTS:block_control>
   $<TARGET_OBJECTS:fft99>
   $<TARGET_OBJECTS:fft>
   $<TARGET_OBJECTS:field_manager>
   $<TARGET_OBJECTS:mosaic>
   $<TARGET_OBJECTS:sat_vapor_pres>
   $<TARGET_OBJECTS:time_manager>
   $<TARGET_OBJECTS:tracer_manager>
   $<TARGET_OBJECTS:astronomy>
   $<TARGET_OBJECTS:column_diagnostics>
   $<TARGET_OBJECTS:diag_manager>
   $<TARGET_OBJECTS:horiz_interp>
   $<TARGET_OBJECTS:oda_tools>
   $<TARGET_OBJECTS:random_numbers>
   $<TARGET_OBJECTS:time_interp>
   $<TARGET_OBJECTS:topography>
   $<TARGET_OBJECTS:amip_interp>
   $<TARGET_OBJECTS:coupler>
   $<TARGET_OBJECTS:data_override>
   $<TARGET_OBJECTS:exchange>
   $<TARGET_OBJECTS:interpolator>
   )
