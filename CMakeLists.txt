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

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# Define the CMake project
project(FMS
  VERSION 2025.03.0
  DESCRIPTION  "GFDL FMS Library"
  HOMEPAGE_URL "https://www.gfdl.noaa.gov/fms"
  LANGUAGES C Fortran)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel|DebugUFS|ReleaseUFS)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo" "DebugUFS" "ReleaseUFS")
endif()

if(NOT CMAKE_C_COMPILER_ID MATCHES "^(Intel|GNU|Clang|IntelLLVM)$")
  message(
    WARNING "Compiler not officially supported: ${CMAKE_C_COMPILER_ID}")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|GNU|IntelLLVM)$")
  message(
    WARNING "Compiler not officially supported: ${CMAKE_Fortran_COMPILER_ID}")
endif()

# Append directory that contains CMake Modules for building FMS
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Build options
option(OPENMP      "Build FMS with OpenMP support"        ON)
option(32BIT       "Build 32-bit (r4) FMS library"        OFF)
option(64BIT       "Build 64-bit (r8) FMS library"        OFF)
option(FPIC        "Build with position independent code" OFF)
option(SHARED_LIBS "Build shared/dynamic libraries"       OFF)

# Options for compiler definitions
option(INTERNAL_FILE_NML     "Enable compiler definition -DINTERNAL_FILE_NML"      ON)
option(ENABLE_QUAD_PRECISION "Enable compiler definition -DENABLE_QUAD_PRECISION"  ON)
option(PORTABLE_KINDS        "Enable compiler definition -DPORTABLE_KINDS"        OFF)
option(GFS_PHYS              "Enable compiler definition -DGFS_PHYS"              OFF)
option(LARGEFILE             "Enable compiler definition -Duse_LARGEFILE"         OFF)
option(WITH_YAML             "Enable compiler definition -Duse_yaml"              OFF)
option(USE_DEPRECATED_IO     "THIS OPTION HAS NO EFFECT AND WILL BE REMOVED IN A FUTURE RELEASE"   OFF)

if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
  set(ENABLE_QUAD_PRECISION OFF)
  message(WARNING "NVHPC does not support building FMS with QUAD_PRECISION!")
endif()

if(32BIT)
  list(APPEND kinds "r4")
endif()
if(64BIT)
  list(APPEND kinds "r8")
endif()
if(NOT kinds)
  message(STATUS "Single Precision 32BIT: ${32BIT}")
  message(STATUS "Double Precision 64BIT: ${64BIT}")
  message(STATUS "No kind specified, building default double precision (with mixed precision support). Real size will not be included in built library name.")
endif()

# Find dependencies
find_package(MPI REQUIRED COMPONENTS C Fortran)
find_package(NetCDF REQUIRED COMPONENTS C Fortran)

# Check for the OpenMP library and set the required compile flags
if (OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS C Fortran)
endif()

if (WITH_YAML)
  find_package(libyaml REQUIRED COMPONENTS C)
  include_directories(${LIBYAML_INCLUDE_DIR})
endif ()

# Enables position independent code (i.e., -fPIC)
if (FPIC)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif ()

# Collect FMS Fortran source files
list(APPEND fms_fortran_src_files
  affinity/fms_affinity.F90
  amip_interp/amip_interp.F90
  astronomy/astronomy.F90
  axis_utils/axis_utils2.F90
  block_control/block_control.F90
  column_diagnostics/column_diagnostics.F90
  constants/constants.F90
  constants/fmsconstants.F90
  constants4/constantsr4.F90
  constants4/fmsconstantsr4.F90
  coupler/atmos_ocean_fluxes.F90
  coupler/gex.F90
  coupler/coupler_types.F90
  coupler/ensemble_manager.F90
  data_override/get_grid_version.F90
  data_override/data_override.F90
  diag_integral/diag_integral.F90
  diag_manager/diag_axis.F90
  diag_manager/diag_data.F90
  diag_manager/diag_grid.F90
  diag_manager/diag_manager.F90
  diag_manager/diag_output.F90
  diag_manager/diag_table.F90
  diag_manager/diag_util.F90
  diag_manager/fms_diag_time_utils.F90
  diag_manager/fms_diag_object.F90
  diag_manager/fms_diag_yaml.F90
  diag_manager/fms_diag_file_object.F90
  diag_manager/fms_diag_field_object.F90
  diag_manager/fms_diag_axis_object.F90
  diag_manager/fms_diag_output_buffer.F90
  diag_manager/fms_diag_input_buffer.F90
  diag_manager/fms_diag_time_reduction.F90
  diag_manager/fms_diag_outfield.F90
  diag_manager/fms_diag_elem_weight_procs.F90
  diag_manager/fms_diag_fieldbuff_update.F90
  diag_manager/fms_diag_bbox.F90
  diag_manager/fms_diag_reduction_methods.F90
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
  field_manager/fm_yaml.F90
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
  monin_obukhov/monin_obukhov_inter.F90
  monin_obukhov/monin_obukhov.F90
  grid_utils/gradient.F90
  mosaic2/grid2.F90
  mosaic2/mosaic2.F90
  mpp/mpp.F90
  mpp/mpp_data.F90
  mpp/mpp_domains.F90
  mpp/mpp_efp.F90
  mpp/mpp_memutils.F90
  mpp/mpp_parameter.F90
  mpp/mpp_utilities.F90
  parser/yaml_parser.F90
  parser/fms_yaml_output.F90
  platform/platform.F90
  random_numbers/mersennetwister.F90
  random_numbers/random_numbers.F90
  sat_vapor_pres/sat_vapor_pres_k.F90
  sat_vapor_pres/sat_vapor_pres.F90
  string_utils/fms_string_utils.F90
  time_interp/time_interp_external2.F90
  time_interp/time_interp.F90
  time_manager/get_cal_time.F90
  time_manager/time_manager.F90
  topography/gaussian_topog.F90
  topography/topography.F90
  tracer_manager/tracer_manager.F90
  tridiagonal/tridiagonal.F90
  libFMS.F90
)

# Collect FMS C source files
list(APPEND fms_c_src_files
  affinity/affinity.c
  fms/fms_stacksize.c
  grid_utils/gradient_c2l.c
  grid_utils/grid_utils.c
  grid_utils/tree_utils.c
  horiz_interp/include/horiz_interp_conserve_xgrid.c
  mpp/mpp_memuse.c
  parser/yaml_parser_binding.c
  parser/yaml_output_functions.c
  string_utils/fms_string_utils_binding.c
  test_fms/fms/test_fmsC.c # TODO just easier to lump this with the rest of the C src
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

# check gettid
include(CheckFunctionExists)
check_function_exists(gettid HAVE_GETTID)
if(HAVE_GETTID)
  list(APPEND fms_defs HAVE_GETTID)
endif()

# Additional (optional) compiler definitions
if(NOT CONSTANTS)
    set(CONSTANTS GFDL)
endif()
if(CONSTANTS STREQUAL "GFS")
  list(APPEND fms_defs GFS_CONSTANTS)
elseif(CONSTANTS STREQUAL "GEOS")
  list(APPEND fms_defs GEOS_CONSTANTS)
elseif(CONSTANTS STREQUAL "GFDL")
  list(APPEND fms_defs GFDL_CONSTANTS)
else()
  message(FATAL_ERROR "CONSTANTS=${CONSTANTS} option not supported")
endif()

if(GFS_PHYS)
  list(APPEND fms_defs GFS_PHYS)
endif()

if(WITH_YAML)
  list(APPEND fms_defs use_yaml)
endif()

if(USE_DEPRECATED_IO)
  message( FATAL_ERROR "FMS_IO HAS BEEN REMOVED! PLEASE UPDATE TO USE FMS2_IO AND REMOVE  "
                   "-DUSE_DEPRECATED_IO=on FROM YOUR OPTIONS OR REVERT TO VERSION 2025.02")
endif()

if(INTERNAL_FILE_NML)
  list(APPEND fms_defs INTERNAL_FILE_NML)
endif()

if(ENABLE_QUAD_PRECISION)
  list(APPEND fms_defs ENABLE_QUAD_PRECISION)
endif()

if(PORTABLE_KINDS)
  list(APPEND fms_defs PORTABLE_KINDS)
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

# If netCDF was not built with HDF5 parallel I/O features, set up the macro -DNO_NC_PARALLEL4
IF(NOT NetCDF_PARALLEL)
  MESSAGE(WARNING "netCDF was not build with HDF5 parallel I/O features, so collective netcdf io is not allowed")
  list(APPEND fms_defs NO_NC_PARALLEL4)
ENDIF()

# TODO remove this section eventually, separate builds are no longer required
foreach(kind ${kinds})

  set(libTgt fms_${kind})
  set(includeDir "include_${kind}")
  set(moduleDir "${CMAKE_CURRENT_BINARY_DIR}/${includeDir}")

  # C
  add_library(${libTgt}_c OBJECT ${fms_c_src_files})

  target_include_directories(${libTgt}_c PRIVATE include
                                                 grid_utils)
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
                                                 fms/include
                                                 fms2_io/include
                                                 string_utils/include
                                                 mpp/include
                                                 column_diagnostics/include
                                                 monin_obukhov/include
                                                 sat_vapor_pres/include
                                                 horiz_interp/include
                                                 diag_integral/include
                                                 random_numbers/include
                                                 diag_manager/include
                                                 constants4
                                                 topography/include
                                                 axis_utils/include
                                                 mosaic2/include
                                                 constants
                                                 astronomy/include
                                                 field_manager/
                                                 field_manager/include
                                                 time_interp/include
                                                 tracer_manager/include
                                                 tridiagonal/include
                                                 interpolator/include
                                                 coupler/include
                                                 data_override/include
                                                 amip_interp/include)
  target_compile_definitions(${libTgt}_f PRIVATE "${fms_defs}")
  target_compile_definitions(${libTgt}_f PRIVATE "${${kind}_defs}")
  set_target_properties(${libTgt}_f PROPERTIES COMPILE_FLAGS ${${kind}_flags})
  set_target_properties(${libTgt}_f PROPERTIES Fortran_MODULE_DIRECTORY
                                               ${moduleDir})
  target_link_libraries(${libTgt}_f PRIVATE NetCDF::NetCDF_Fortran
                                            MPI::MPI_Fortran)
  if(OpenMP_Fortran_FOUND)
    target_link_libraries(${libTgt}_f PRIVATE OpenMP::OpenMP_Fortran)
  endif()


  # Check if gnu 10 or higher
  # this should only be needed with mpich, but wasn't able to find a good way to find the MPI flavor consistently
  if ( CMAKE_Fortran_COMPILER_VERSION MATCHES "1[0-9]\.[0-9]*\.[0-9]*" AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    include(CheckFortranCompilerFlag)
    check_fortran_compiler_flag("-fallow-argument-mismatch" _arg_mismatch_flag)
    if(_arg_mismatch_flag)
      message(STATUS "Adding -fallow-argument-mismatch flag to compile with GCC >=10 and MPICH")
      target_compile_options(${libTgt}_f PRIVATE "-fallow-argument-mismatch;-w")
    endif()
  endif()

  if ( CMAKE_Fortran_COMPILER_VERSION MATCHES "2025\.1\..*" AND CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    message(FATAL_ERROR "Intel LLVM Fortran compiler version 2025.1 is not supported due to internal compiler error")
  endif()

  # FMS (C + Fortran)
  if (SHARED_LIBS)
      message(STATUS "Shared library target: ${libTgt}")
      add_library(${libTgt} SHARED $<TARGET_OBJECTS:${libTgt}_c>
                                   $<TARGET_OBJECTS:${libTgt}_f>)
  else ()
      message(STATUS "Static library target: ${libTgt}")
      add_library(${libTgt} STATIC $<TARGET_OBJECTS:${libTgt}_c>
                                   $<TARGET_OBJECTS:${libTgt}_f>)
  endif ()

  target_include_directories(${libTgt} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fms>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fms/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fms2_io/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/column_diagnostics/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/astronomy/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/monin_obukhov/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sat_vapor_pres/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/field_manager/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/grid_utils>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/horiz_interp/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mosaic2/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/string_utils/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpp/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/diag_manager/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/topography/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/random_numbers/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/diag_integral/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/diag_manager/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tracer_manager/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interpolator/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/coupler/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mosaic2/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/data_override/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tridiagonal/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/amip_interp/include>)
  set(moduleDir "${CMAKE_CURRENT_BINARY_DIR}/${includeDir}")
  target_include_directories(${libTgt} INTERFACE
    $<BUILD_INTERFACE:${moduleDir}>
    $<INSTALL_INTERFACE:${includeDir}>)

  target_compile_definitions(${libTgt} PRIVATE "${fms_defs}")
  target_compile_definitions(${libTgt} PRIVATE "${${kind}_defs}")

  target_link_libraries(${libTgt} PUBLIC NetCDF::NetCDF_C
                                         NetCDF::NetCDF_Fortran
                                         MPI::MPI_Fortran)

  if(OpenMP_Fortran_FOUND)
    target_link_libraries(${libTgt} PRIVATE OpenMP::OpenMP_C OpenMP::OpenMP_Fortran)
  endif()

  add_library(FMS::${libTgt} ALIAS ${libTgt})

  list(APPEND LIB_TARGETS ${libTgt})
  install(DIRECTORY ${moduleDir}    DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(FILES ${fms_header_files} DESTINATION ${CMAKE_INSTALL_PREFIX}/${includeDir})

endforeach()

# default build, same as an r8 build just different name
if(NOT kinds)
    message(STATUS "building default kind")
    set(libTgt "fms")
    set(includeDir "include")
    set(moduleDir "${CMAKE_CURRENT_BINARY_DIR}/${includeDir}")
    set(kind "r8")

    # C
    add_library(${libTgt}_c OBJECT ${fms_c_src_files})

    target_include_directories(${libTgt}_c PRIVATE include
                                                   grid_utils)
    target_compile_definitions(${libTgt}_c PRIVATE "${fms_defs}")
    target_link_libraries(${libTgt}_c PRIVATE NetCDF::NetCDF_C
                                              MPI::MPI_C)

    if(OpenMP_C_FOUND)
      target_link_libraries(${libTgt}_c PRIVATE OpenMP::OpenMP_C)
    endif()

    if(WITH_YAML)
      target_link_libraries(${libTgt}_c PRIVATE libyaml_C)
      set_target_properties(${libTgt}_c PROPERTIES COMPILE_FLAGS "-L${LIBYAML_LIBRARIES} -lyaml")
    endif()


    # Fortran
    add_library(${libTgt}_f OBJECT ${fms_fortran_src_files})

    target_include_directories(${libTgt}_f PRIVATE include
                                                   fms
                                                   fms/include
                                                   fms2_io/include
                                                   string_utils/include
                                                   mpp/include
                                                   column_diagnostics/include
                                                   monin_obukhov/include
                                                   sat_vapor_pres/include
                                                   horiz_interp/include
                                                   diag_integral/include
                                                   random_numbers/include
                                                   diag_manager/include
                                                   constants4
                                                   topography/include
                                                   axis_utils/include
                                                   mosaic2/include
                                                   constants
                                                   astronomy/include
                                                   field_manager/include
                                                   time_interp/include
                                                   tracer_manager/include
                                                   tridiagonal/include
                                                   interpolator/include
                                                   coupler/include
                                                   data_override/include
                                                   amip_interp/include)
    target_compile_definitions(${libTgt}_f PRIVATE "${fms_defs}")
    target_compile_definitions(${libTgt}_f PRIVATE "${${kind}_defs}")
    set_target_properties(${libTgt}_f PROPERTIES COMPILE_FLAGS ${${kind}_flags})
    set_target_properties(${libTgt}_f PROPERTIES Fortran_MODULE_DIRECTORY
                                                 ${moduleDir})
    target_link_libraries(${libTgt}_f PRIVATE NetCDF::NetCDF_Fortran
                                              MPI::MPI_Fortran)
    if(OpenMP_Fortran_FOUND)
      target_link_libraries(${libTgt}_f PRIVATE OpenMP::OpenMP_Fortran)
    endif()

    # Check if gnu 10 or higher
    # this should only be needed with mpich, but wasn't able to find a good way to find the MPI flavor consistently
    if ( CMAKE_Fortran_COMPILER_VERSION MATCHES "1[0-9]\.[0-9]*\.[0-9]*" AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      include(CheckFortranCompilerFlag)
      check_fortran_compiler_flag("-fallow-argument-mismatch" _arg_mismatch_flag)
      if(_arg_mismatch_flag)
        message(STATUS "Adding -fallow-argument-mismatch flag to compile with GCC >=10 and MPICH")
        target_compile_options(${libTgt}_f PRIVATE "-fallow-argument-mismatch;-w")
      endif()
    endif()

    # FMS (C + Fortran)
    if (SHARED_LIBS)
        message(STATUS "Shared library target: ${libTgt}")
        add_library(${libTgt} SHARED $<TARGET_OBJECTS:${libTgt}_c>
                                     $<TARGET_OBJECTS:${libTgt}_f>)
    else ()
        message(STATUS "Static library target: ${libTgt}")
        add_library(${libTgt} STATIC $<TARGET_OBJECTS:${libTgt}_c>
                                     $<TARGET_OBJECTS:${libTgt}_f>)
    endif ()

    target_include_directories(${libTgt} PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fms>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fms/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fms2_io/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/column_diagnostics/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/astronomy/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/monin_obukhov/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sat_vapor_pres/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/field_manager/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/grid_utils>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/horiz_interp/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mosaic2/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/string_utils/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mpp/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/diag_manager/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/topography/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/random_numbers/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/diag_integral/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/diag_manager/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tracer_manager/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interpolator/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/coupler/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mosaic2/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/data_override/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tridiagonal/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/amip_interp/include>)

    set(moduleDir "${CMAKE_CURRENT_BINARY_DIR}/${includeDir}")
    target_include_directories(${libTgt} INTERFACE
      $<BUILD_INTERFACE:${moduleDir}>
      $<INSTALL_INTERFACE:${includeDir}>)

    target_compile_definitions(${libTgt} PRIVATE "${fms_defs}")
    target_compile_definitions(${libTgt} PRIVATE "${${kind}_defs}")
    set_target_properties(${libTgt}_f PROPERTIES COMPILE_FLAGS ${${kind}_flags})

    target_link_libraries(${libTgt} PUBLIC NetCDF::NetCDF_C
                                           NetCDF::NetCDF_Fortran
                                           MPI::MPI_Fortran)

    if(OpenMP_Fortran_FOUND)
      target_link_libraries(${libTgt} PRIVATE OpenMP::OpenMP_C OpenMP::OpenMP_Fortran)
    endif()

    add_library(FMS::${libTgt} ALIAS ${libTgt})

    list(APPEND LIB_TARGETS ${libTgt})
    install(DIRECTORY ${moduleDir}    DESTINATION ${CMAKE_INSTALL_PREFIX})
    install(FILES ${fms_header_files} DESTINATION ${CMAKE_INSTALL_PREFIX}/${includeDir})


endif()

install(
  TARGETS ${LIB_TARGETS}
  EXPORT FMSExports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

### Unit Testing
include(CTest)

# TODO autotools also checks if srun is available
set(MPI_LAUNCHER "mpirun")
# used in the test-lib.sh.in to make it behave differently when parsed by cmake
set(USING_CMAKE "true")

# parse and add build info to test script util file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_fms/test-lib.sh.in ${CMAKE_CURRENT_SOURCE_DIR}/test-lib.sh
               @ONLY)

# build any helper modules that are used by unit tests
list(APPEND TEST_MODS_SRC
  test_fms/diag_manager/testing_utils.F90
  test_fms/fms/test_fms.F90
  test_fms/fms2_io/argparse.F90
  test_fms/fms2_io/setup.F90
  test_fms/mosaic2/write_files.F90
  test_fms/mpp/compare_data_checksums.F90
  test_fms/mpp/compare_data_checksums_int.F90
  test_fms/mpp/fill_halo.F90
  test_fms/mpp/test_domains_utility_mod.F90
  test_fms/mpp/test_mpp_update_domains_int.F90
  test_fms/mpp/test_mpp_update_domains_real.F90
  test_fms/mpp/test_mpp_update_domains_ad.F90
  test_fms/mpp/test_system_clock.F90)

add_library(testLibs "${TEST_MODS_SRC}")
target_link_libraries(testLibs FMS::fms)
target_compile_definitions(testLibs PRIVATE TEST_MOS_KIND_=r8_kind)
target_include_directories(testLibs PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include
                                    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_fms/fms/include
                                    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_fms/data_override/include
                                    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_fms/mosaic2
                                    PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(testLibs PROPERTIES Fortran_MODULE_DIRECTORY
                                               ${moduleDir})

# add all our test source code

# this lists all source code that creates a single executable for unit testing (r8 default)
list(APPEND TEST_SRC_SINGLE_TARGET
  test_fms/affinity/test_affinity.F90
  test_fms/block_control/test_block_control.F90
  test_fms/coupler/test_gex.F90
  test_fms/data_override/test_data_override_init.F90
  test_fms/data_override/test_data_override_ongrid.F90
  test_fms/diag_manager/test_output_every_freq.F90
  test_fms/diag_manager/test_dm_weights.F90
  test_fms/diag_manager/test_diag_manager.F90
  test_fms/diag_manager/test_diag_manager_time.F90
  test_fms/diag_manager/test_diag_update_buffer.F90
  test_fms/diag_manager/test_diag_yaml.F90
  test_fms/diag_manager/test_diag_ocean.F90
  test_fms/diag_manager/test_modern_diag.F90
  test_fms/diag_manager/test_diag_buffer.F90
  test_fms/diag_manager/test_flexible_time.F90
  test_fms/diag_manager/test_flush_nc_file.F90
  test_fms/diag_manager/test_diag_out_yaml.F90
  test_fms/diag_manager/test_reduction_methods.F90
  test_fms/diag_manager/testing_utils.F90
  test_fms/diag_manager/test_diag_diurnal.F90
  test_fms/diag_manager/check_time_none.F90
  test_fms/diag_manager/check_time_min.F90
  test_fms/diag_manager/check_time_max.F90
  test_fms/diag_manager/check_time_avg.F90
  test_fms/diag_manager/check_time_diurnal.F90
  test_fms/diag_manager/check_time_pow.F90
  test_fms/diag_manager/check_time_rms.F90
  test_fms/diag_manager/test_cell_measures.F90
  test_fms/diag_manager/check_subregional.F90
  test_fms/diag_manager/test_var_masks.F90
  test_fms/diag_manager/check_var_masks.F90
  test_fms/diag_manager/test_multiple_send_data.F90
  test_fms/diag_manager/test_prepend_date.F90
  test_fms/diag_manager/test_diag_multi_file.F90
  test_fms/diag_manager/test_ens_runs.F90
  test_fms/diag_manager/test_diag_attribute_add.F90
  test_fms/diag_manager/check_new_file_freq.F90
  test_fms/drifters/test_cloud_interpolator.F90
  test_fms/drifters/test_drifters_io.F90
  test_fms/drifters/test_drifters_input.F90
  test_fms/drifters/test_drifters_comm.F90
  test_fms/drifters/test_drifters_core.F90
  test_fms/drifters/test_quicksort.F90
  test_fms/drifters/test_drifters.F90
  test_fms/exchange/test_xgrid.F90
  test_fms/field_manager/test_field_table_read.F90
  test_fms/fms2_io/test_get_is_valid.F90
  test_fms/fms2_io/test_fms2_io.F90
  test_fms/fms2_io/test_atmosphere_io.F90
  test_fms/fms2_io/test_io_simple.F90
  test_fms/fms2_io/setup.F90
  test_fms/fms2_io/argparse.F90
  test_fms/fms2_io/test_global_att.F90
  test_fms/fms2_io/test_io_with_mask.F90
  test_fms/fms2_io/test_bc_restart.F90
  test_fms/fms2_io/test_get_mosaic_tile_grid.F90
  test_fms/fms2_io/test_read_ascii_file.F90
  test_fms/fms2_io/test_file_appendix.F90
  test_fms/fms2_io/test_unlimit_compressed.F90
  test_fms/fms2_io/test_chunksizes.F90
  test_fms/fms2_io/test_packed_reads.F90
  test_fms/fms2_io/test_compressed_writes.F90
  test_fms/fms2_io/test_domain_io.F90
  test_fms/fms2_io/test_collective_io.F90
  test_fms/fms2_io/test_generalized_indices.F90
  test_fms/fms/test_fms.F90
  test_fms/interpolator/test_interpolator.F90
  test_fms/mpp/test_clock_init.F90
  test_fms/mpp/test_domains_simple.F90
  test_fms/mpp/test_domains_utility_mod.F90
  test_fms/mpp/test_global_arrays.F90
  test_fms/mpp/test_minmax.F90
  test_fms/mpp/test_mpp_alltoall.F90
  test_fms/mpp/test_mpp_broadcast.F90
  test_fms/mpp/test_mpp_chksum.F90
  test_fms/mpp/test_mpp_clock_begin_end_id.F90
  test_fms/mpp/test_mpp_domains.F90
  test_fms/mpp/test_mpp.F90
  test_fms/mpp/test_mpp_gatscat.F90
  test_fms/mpp/test_mpp_get_ascii_lines.F90
  test_fms/mpp/test_mpp_global_field.F90
  test_fms/mpp/test_mpp_global_field_ug.F90
  test_fms/mpp/test_mpp_global_sum_ad.F90
  test_fms/mpp/test_mpp_init_logfile.F90
  test_fms/mpp/test_mpp_mem_dump.F90
  test_fms/mpp/test_mpp_memuse.F90
  test_fms/mpp/test_mpp_memutils_begin_2x.F90
  test_fms/mpp/test_mpp_memutils_begin_end.F90
  test_fms/mpp/test_mpp_memutils_end_before_begin.F90
  test_fms/mpp/test_mpp_nesting.F90
  test_fms/mpp/test_mpp_npes.F90
  test_fms/mpp/test_mpp_pe.F90
  test_fms/mpp/test_mpp_print_memuse_stats_file.F90
  test_fms/mpp/test_mpp_print_memuse_stats_stderr.F90
  test_fms/mpp/test_mpp_root_pe.F90
  test_fms/mpp/test_mpp_sendrecv.F90
  test_fms/mpp/test_mpp_sum.F90
  test_fms/mpp/test_mpp_transmit.F90
  test_fms/mpp/test_mpp_update_domains.F90
  test_fms/mpp/test_mpp_update_domains_ad.F90
  test_fms/mpp/test_peset.F90
  test_fms/mpp/test_mpp_read_ascii_file.F90
  test_fms/mpp/test_read_input_nml.F90
  test_fms/mpp/test_redistribute_int.F90
  test_fms/mpp/test_stderr.F90
  test_fms/mpp/test_stdin.F90
  test_fms/mpp/test_stdlog.F90
  test_fms/mpp/test_stdout.F90
  test_fms/mpp/test_super_grid.F90
  test_fms/mpp/test_system_clock.F90
  test_fms/mpp/test_update_domains_performance.F90
  test_fms/parser/test_yaml_parser.F90
  test_fms/parser/check_crashes.F90
  test_fms/parser/parser_demo.F90
  test_fms/parser/parser_demo2.F90
  test_fms/parser/test_output_yaml.F90
  test_fms/parser/generic_blocks.F90
  test_fms/string_utils/test_string_utils.F90
  test_fms/time_manager/test_time_manager.F90
  )
# list any test sources that are used for two executables (r4 and r8)
list(APPEND TEST_SRC_MULTI_TARGET
  test_fms/astronomy/test_astronomy.F90
  test_fms/axis_utils/test_axis_utils.F90
  test_fms/column_diagnostics/test_column_diagnostics.F90
  test_fms/coupler/test_coupler_2d.F90
  test_fms/coupler/test_coupler_3d.F90
  test_fms/coupler/test_coupler_types.F90
  test_fms/coupler/test_atmos_ocean_fluxes.F90
  test_fms/data_override/test_data_override.F90
  test_fms/data_override/test_get_grid_v1.F90
  test_fms/diag_integral/test_diag_integral.F90
  test_fms/field_manager/test_field_manager.F90
  test_fms/horiz_interp/test_horiz_interp.F90
  test_fms/horiz_interp/test_create_xgrid_order2.F90
  test_fms/interpolator/test_interpolator2.F90
  test_fms/monin_obukhov/test_monin_obukhov.F90
  test_fms/mosaic2/test_mosaic2.F90
  test_fms/mosaic2/test_grid2.F90
  test_fms/random_numbers/test_random_numbers.F90
  test_fms/sat_vapor_pres/test_sat_vapor_pres.F90
  test_fms/time_interp/test_time_interp.F90
  test_fms/time_interp/test_time_interp_external.F90
  test_fms/topography/test_topography.F90
  test_fms/tracer_manager/test_tracer_manager.F90
  test_fms/tridiagonal/test_tridiagonal.F90
  )

# add all the test executables and set up linking + compile flags
foreach (testFile ${TEST_SRC_SINGLE_TARGET})
  get_filename_component (TName ${testFile} NAME_WE)

  add_executable(${TName} ${testFile})
  target_compile_definitions(${TName} PRIVATE "${r8_defs}")
  set_target_properties(${TName} PROPERTIES COMPILE_FLAGS ${r8_flags})
  target_link_libraries(${TName} PUBLIC FMS::fms
                                 PRIVATE testLibs
                                   )

  if(WITH_YAML)
    target_link_libraries(${TName} PRIVATE yaml ${LIBYAML_LIBRARIES})
  endif()

  if(OPENMP)
    target_link_libraries(${TName} PRIVATE OpenMP::OpenMP_Fortran)
  endif()

  set_target_properties(${TName} PROPERTIES Fortran_MODULE_DIRECTORY
                                            ${moduleDir})
  target_include_directories(${TName} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include
                                       PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_fms/fms/include
                                       PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_fms/data_override/include
                                       PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_fms/mosaic2
                                       PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
endforeach ()

foreach (testFile ${TEST_SRC_MULTI_TARGET})
  get_filename_component (TName ${testFile} NAME_WE)

  add_executable(${TName}_r8 ${testFile})
  target_compile_definitions(${TName}_r8 PRIVATE "${r8_defs}")
  set_target_properties(${TName}_r8 PROPERTIES COMPILE_FLAGS ${r8_flags})
  target_link_libraries(${TName}_r8 PUBLIC FMS::fms
                                    PRIVATE testLibs
                                    )

  add_executable(${TName}_r4 ${testFile})
  target_compile_definitions(${TName}_r4 PRIVATE "${r4_defs}")
  # seems counterintuitive but r4 tests use r8 default
  # they specify kind values explicitly with the preprocessor where needed (TEST_FMS_KIND_)
  set_target_properties(${TName}_r4 PROPERTIES COMPILE_FLAGS ${r8_flags})
  target_link_libraries(${TName}_r4 PUBLIC FMS::fms
                                    PRIVATE testLibs
                                    )

  if(WITH_YAML)
    target_link_libraries(${TName}_r4 PRIVATE yaml ${LIBYAML_LIBRARIES})
    target_link_libraries(${TName}_r8 PRIVATE yaml ${LIBYAML_LIBRARIES})
  endif()

  if(OPENMP)
    target_link_libraries(${TName}_r4 PRIVATE OpenMP::OpenMP_Fortran)
    target_link_libraries(${TName}_r8 PRIVATE OpenMP::OpenMP_Fortran)
  endif()

  set_target_properties(${TName}_r4 PROPERTIES Fortran_MODULE_DIRECTORY
                                            ${moduleDir})
  set_target_properties(${TName}_r8 PROPERTIES Fortran_MODULE_DIRECTORY
                                            ${moduleDir})
  target_include_directories(${TName}_r4 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include
                                       PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_fms/fms/include
                                       PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_fms/data_override/include
                                       PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_fms/mosaic2
                                       PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
  target_include_directories(${TName}_r8 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include
                                       PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_fms/fms/include
                                       PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_fms/data_override/include
                                       PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test_fms/mosaic2
                                       PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
  #set macros used to test mixed precision
  target_compile_definitions(${TName}_r8 PRIVATE TEST_FMS_KIND_=r8_kind)
  target_compile_definitions(${TName}_r4 PRIVATE TEST_FMS_KIND_=r4_kind)
  target_compile_definitions(${TName}_r8 PRIVATE FMS_TEST_BC_TYPE_=bc)
  target_compile_definitions(${TName}_r4 PRIVATE FMS_TEST_BC_TYPE_=bc_r4)
  target_compile_definitions(${TName}_r8 PRIVATE FMS_TEST_BC_TYPE_=bc)
  target_compile_definitions(${TName}_r4 PRIVATE FMS_TEST_BC_TYPE_=bc_r4)
  target_compile_definitions(${TName}_r8 PRIVATE WRITE_FILES_MOD_=write_files_r8)
  target_compile_definitions(${TName}_r4 PRIVATE WRITE_FILES_MOD_=write_files_r4)
endforeach ()

# gather all test scripts
file(GLOB_RECURSE TEST_SCRIPTS ${CMAKE_CURRENT_SOURCE_DIR}/test_fms/*/*.sh)
# leave out scripts that aren't tests
list(FILTER TEST_SCRIPTS EXCLUDE REGEX "check_crashes.*")
list(FILTER TEST_SCRIPTS EXCLUDE REGEX "test-lib*")
list(FILTER TEST_SCRIPTS EXCLUDE REGEX "tap-driver.*")

# add scripts to ctest
foreach (testScript ${TEST_SCRIPTS})
  get_filename_component(testName ${testScript} NAME_WE)
  get_filename_component(testDirNameFull ${testScript} DIRECTORY)
  get_filename_component(testDirName ${testDirNameFull} NAME_WE)
  add_test(NAME ${testName}
           COMMAND ${testScript}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  # all tests are labeled with directory name
  set_property(TEST ${testName} PROPERTY LABELS ${testDirName})
  # set an environment variable so the scripts know to skip tests that use libyaml
  if(NOT WITH_YAML)
    set_tests_properties(${testName} PROPERTIES
        ENVIRONMENT "parser_skip=skip"
    )
  endif()
endforeach()

set(CMAKE_CTEST_ARGUMENTS "--output-on-failure")

### Package config
include(CMakePackageConfigHelpers)
set(CONFIG_INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fms)

export(EXPORT FMSExports
  NAMESPACE FMS::
  FILE fms-targets.cmake)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FMSConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/fms-config.cmake
  INSTALL_DESTINATION ${CONFIG_INSTALL_DESTINATION})
install(FILES ${CMAKE_SOURCE_DIR}/cmake/FindNetCDF.cmake ${CMAKE_CURRENT_BINARY_DIR}/fms-config.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/fms-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/fms-config-version.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

install(EXPORT FMSExports
  NAMESPACE FMS::
  FILE fms-targets.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

# pkgconf
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix ${CMAKE_INSTALL_PREFIX})
set(libdir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(includedir ${CMAKE_INSTALL_PREFIX}/${includeDir})

set(CC ${CMAKE_C_COMPILER})
set(FC ${CMAKE_Fortran_COMPILER})
set(CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}")
set(CPPFLAGS "${CMAKE_CPP_FLAGS} ${CMAKE_CPP_FLAGS_${CMAKE_BUILD_TYPE}}")
set(FCFLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}")
set(LDFLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_${CMAKE_BUILD_TYPE}}")

set(VERSION ${PROJECT_VERSION})

# TODO: If FMS depends on a library that is built as a static library, it
#       should be listed here as an ldflag.
set(LIBS "")

if(NOT ${NetCDF_Fortran_LIBRARY_SHARED})
  # autotools: Libs.private: -lnetcdff -lnetcdf
  string(APPEND LIBS ${NetCDF_Fortran_LIBRARIES})
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/FMS.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/FMS.pc @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/FMS.pc
              DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
              COMPONENT utilities)
