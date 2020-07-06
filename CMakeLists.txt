#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# Copyright (c) GFDL, @underwoo

cmake_minimum_required(VERSION 3.7...3.15 FATAL_ERROR)

# Define the CMake project
project(FMS
  VERSION 2020.02
  DESCRIPTION  "GFDL FMS Library"
  HOMEPAGE_URL "https://www.gfdl.noaa.gov/fms"
  LANGUAGES C Fortran)

include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT CMAKE_C_COMPILER_ID MATCHES "^(Intel|GNU|Clang)$")
  message(
    WARNING "Compiler not officially supported: ${CMAKE_C_COMPILER_ID}")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|GNU)$")
  message(
    WARNING "Compiler not officially supported: ${CMAKE_Fortran_COMPILER_ID}")
endif()

# Append directory that contains CMake Modules for building FMS
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Build options
option(OPENMP    "Build FMS with OpenMP support" OFF)
option(32BIT     "Build 32-bit (r4) FMS library"  ON)
option(64BIT     "Build 64-bit (r8) FMS library" OFF)

# Options for compiler definitions
option(INTERNAL_FILE_NML     "Enable compiler definition -DINTERNAL_FILE_NML"      ON)
option(ENABLE_QUAD_PRECISION "Enable compiler definition -DENABLE_QUAD_PRECISION"  ON)
option(GFS_PHYS              "Enable compiler definition -DGFS_PHYS"              OFF)
option(LARGEFILE             "Enable compiler definition -Duse_LARGEFILE"         OFF)

if(32BIT)
  list(APPEND kinds "r4")
endif()
if(64BIT)
  list(APPEND kinds "r8")
endif()
if(NOT kinds)
  message(STATUS "Single Precision 32BIT: ${32BIT}")
  message(STATUS "Double Precision 64BIT: ${64BIT}")
  message(FATAL_ERROR "Either 32BIT or 64BIT should be ON")
endif()

# Find dependencies
find_package(MPI REQUIRED COMPONENTS C Fortran)
find_package(NetCDF REQUIRED COMPONENTS C Fortran)

# Check for the OpenMP library and set the required compile flags
if (OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS C Fortran)
endif()

# Collect FMS Fortran source files
list(APPEND fms_fortran_src_files
  affinity/fms_affinity.F90
  amip_interp/amip_interp.F90
  astronomy/astronomy.F90
  axis_utils/axis_utils.F90
  axis_utils/axis_utils2.F90
  block_control/block_control.F90
  column_diagnostics/column_diagnostics.F90
  constants/constants.F90
  coupler/atmos_ocean_fluxes.F90
  coupler/coupler_types.F90
  coupler/ensemble_manager.F90
  data_override/data_override.F90
  diag_integral/diag_integral.F90
  diag_manager/diag_axis.F90
  diag_manager/diag_data.F90
  diag_manager/diag_grid.F90
  diag_manager/diag_manager.F90
  diag_manager/diag_output.F90
  diag_manager/diag_table.F90
  diag_manager/diag_util.F90
  drifters/cloud_interpolator.F90
  drifters/drifters.F90
  drifters/drifters_comm.F90
  drifters/drifters_core.F90
  drifters/drifters_input.F90
  drifters/drifters_io.F90
  drifters/quicksort.F90
  exchange/stock_constants.F90
  exchange/xgrid.F90
  field_manager/field_manager.F90
  field_manager/fm_util.F90
  fms/fms_io.F90
  fms/fms.F90
  fms2_io/blackboxio.F90
  fms2_io/fms_io_utils.F90
  fms2_io/fms_netcdf_domain_io.F90
  fms2_io/fms_netcdf_unstructured_domain_io.F90
  fms2_io/fms2_io.F90
  fms2_io/netcdf_io.F90
  horiz_interp/horiz_interp_bicubic.F90
  horiz_interp/horiz_interp_bilinear.F90
  horiz_interp/horiz_interp_conserve.F90
  horiz_interp/horiz_interp_spherical.F90
  horiz_interp/horiz_interp_type.F90
  horiz_interp/horiz_interp.F90
  interpolator/interpolator.F90
  memutils/memutils.F90
  monin_obukhov/monin_obukhov_kernel.F90
  monin_obukhov/monin_obukhov.F90
  mosaic/gradient.F90
  mosaic/grid.F90
  mosaic/mosaic.F90
  mosaic/mosaic2.F90
  mpp/mpp.F90
  mpp/mpp_data.F90
  mpp/mpp_domains.F90
  mpp/mpp_efp.F90
  mpp/mpp_io.F90
  mpp/mpp_memutils.F90
  mpp/mpp_parameter.F90
  mpp/mpp_pset.F90
  mpp/mpp_utilities.F90
  platform/platform.F90
  random_numbers/MersenneTwister.F90
  random_numbers/random_numbers.F90
  sat_vapor_pres/sat_vapor_pres_k.F90
  sat_vapor_pres/sat_vapor_pres.F90
  time_interp/time_interp_external.F90
  time_interp/time_interp_external2.F90
  time_interp/time_interp.F90
  time_manager/get_cal_time.F90
  time_manager/time_manager.F90
  topography/gaussian_topog.F90
  topography/topography.F90
  tracer_manager/tracer_manager.F90
  tridiagonal/tridiagonal.F90
)

# Collect FMS C source files
list(APPEND fms_c_src_files
  affinity/affinity.c
  mosaic/create_xgrid.c
  mosaic/gradient_c2l.c
  mosaic/interp.c
  mosaic/mosaic_util.c
  mosaic/read_mosaic.c
  mpp/mpp_memuse.c
  mpp/nsclock.c
  mpp/threadloc.c
)

# Collect FMS header files
list(APPEND fms_header_files
  include/file_version.h
  include/fms_platform.h
)

# Standard FMS compiler definitions
list(APPEND fms_defs
  use_libMPI
  use_netCDF)

# Additional (optional) compiler definitions
if(GFS_PHYS)
  list(APPEND fms_defs GFS_PHYS)
endif()

if(INTERNAL_FILE_NML)
  list(APPEND fms_defs INTERNAL_FILE_NML)
endif()

if(ENABLE_QUAD_PRECISION)
  list(APPEND fms_defs ENABLE_QUAD_PRECISION)
endif()

if(LARGEFILE)
  list(APPEND fms_defs use_LARGEFILE)
endif()

# Precision-based compiler definitions
if(32BIT)
  list(APPEND r4_defs OVERLOAD_R4 OVERLOAD_R8)
endif()

# Add platform specific compiler definitions
if(APPLE)
  list(APPEND fms_defs __APPLE__)
endif()

# Obtain compiler-specific flags
include(fms_compiler_flags)

foreach(kind ${kinds})

  set(libTgt fms_${kind})
  set(includeDir "include_${kind}")
  set(moduleDir "${CMAKE_CURRENT_BINARY_DIR}/${includeDir}")

  # C
  add_library(${libTgt}_c OBJECT ${fms_c_src_files})

  target_include_directories(${libTgt}_c PRIVATE include)
  target_compile_definitions(${libTgt}_c PRIVATE "${fms_defs}")

  target_link_libraries(${libTgt}_c PRIVATE NetCDF::NetCDF_C
                                            MPI::MPI_C)

  if(OpenMP_C_FOUND)
    target_link_libraries(${libTgt}_c PRIVATE OpenMP::OpenMP_C)
  endif()

  # Fortran
  add_library(${libTgt}_f OBJECT ${fms_fortran_src_files})

  target_include_directories(${libTgt}_f PRIVATE include
                                                 fms
                                                 fms2_io/include
                                                 mpp/include)
  target_compile_definitions(${libTgt}_f PRIVATE "${fms_defs}")
  target_compile_definitions(${libTgt}_f PRIVATE "${${kind}_defs}")

  set_target_properties(${libTgt}_f PROPERTIES COMPILE_FLAGS
                                               "${${kind}_flags}")
  set_target_properties(${libTgt}_f PROPERTIES Fortran_MODULE_DIRECTORY
                                               ${moduleDir})

  target_link_libraries(${libTgt}_f PRIVATE NetCDF::NetCDF_Fortran
                                            MPI::MPI_Fortran)

  if(OpenMP_Fortran_FOUND)
    target_link_libraries(${libTgt}_f PRIVATE OpenMP::OpenMP_Fortran)
  endif()

  # FMS (C + Fortran)
  add_library(${libTgt} STATIC $<TARGET_OBJECTS:${libTgt}_c>
                               $<TARGET_OBJECTS:${libTgt}_f>)

  target_include_directories(${libTgt} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fms>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fms2_io/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpp/include>)

  target_include_directories(${libTgt} INTERFACE
    $<BUILD_INTERFACE:${moduleDir}>
    $<INSTALL_INTERFACE:${includeDir}>)

  target_compile_definitions(${libTgt} PUBLIC "${fms_defs}")
  target_compile_definitions(${libTgt} PUBLIC "${${kind}_defs}")

  target_link_libraries(${libTgt} PUBLIC NetCDF::NetCDF_Fortran
                                         MPI::MPI_Fortran)

  if(OpenMP_Fortran_FOUND)
    target_link_libraries(${libTgt} PUBLIC OpenMP::OpenMP_Fortran)
  endif()

  add_library(FMS::${libTgt} ALIAS ${libTgt})

  list(APPEND LIB_TARGETS ${libTgt})
  install(DIRECTORY ${moduleDir}    DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(FILES ${fms_header_files} DESTINATION ${CMAKE_INSTALL_PREFIX}/${includeDir})

endforeach()

install(
  TARGETS ${LIB_TARGETS}
  EXPORT FMSExports
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

### Package config
include(CMakePackageConfigHelpers)
set(CONFIG_INSTALL_DESTINATION lib/cmake/fms)

export(EXPORT FMSExports
  NAMESPACE FMS::
  FILE fms-targets.cmake)

configure_package_config_file(
  ${CMAKE_SOURCE_DIR}/cmake/FMSConfig.cmake.in ${CMAKE_BINARY_DIR}/fms-config.cmake
  INSTALL_DESTINATION ${CONFIG_INSTALL_DESTINATION})
install(FILES ${CMAKE_BINARY_DIR}/fms-config.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

write_basic_package_version_file(
  ${CMAKE_BINARY_DIR}/fms-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_BINARY_DIR}/fms-config-version.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

install(EXPORT FMSExports
  NAMESPACE FMS::
  FILE fms-targets.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})
