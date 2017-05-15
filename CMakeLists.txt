# -*-cmake-*-

project(libfms C Fortran)

set(CMAKE_MODULE_PATH ${libfms_SOURCE_DIR}/CMakeModules ${CMAKE_MODULE_PATH})

# Minimul version of CMake required
cmake_minimum_required (VERSION 3.7)

# Compiler Flags
set(GNU_F90_FLAGS_BASE "-fcray-pointer -fdefault-double-8 -fdefault-real-8 -Waliasing -ffree-line-length-none -fno-range-check")
set(GNU_F90_CPPFLAGS "-DINTERNAL_FILE_NML")

set(CMAKE_Fortran_FLAGS "${GNU_F90_FLAGS_BASE}")

# Required external packages/libraries
find_package(NetCDF)
find_package(NetCDFF REQUIRED)
find_package(MPI)

# Include directory required by all
include_directories("${PROJECT_SOURCE_DIR}/include")

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

add_library(libfms
   $<TARGET_OBJECTS:constants>
   $<TARGET_OBJECTS:platform>
   $<TARGET_OBJECTS:tridiagonal>
   $<TARGET_OBJECTS:mpp>
   $<TARGET_OBJECTS:memutils>
   $<TARGET_OBJECTS:drifters>
   $<TARGET_OBJECTS:fms>
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

if (NETCDF_FOUND OR NETCDFF_FOUND)
  target_compile_definitions(libfms PUBLIC "use_netCDF")
endif (NETCDF_FOUND OR NETCDFF_FOUND)

if (NETCDF_FOUND)
  include_directories(${NETCDF_INCLUDE_PATH})
endif (NETCDF_FOUND)
if (NETCDFF_FOUND)
  include_directories(${NETCDFF_INCLUDE_PATH})
endif (NETCDFF_FOUND)

if (MPI_C_FOUND OR MPI_Fortran_FOUND)
  target_compile_definitions(libfms PUBLIC "use_MIP")
endif (MPI_C_FOUND OR MPI_Fortran_FOUND)

target_compile_definitions(libfms PUBLIC "INTERNAL_FILE_NML")
