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

# **********************************
# amip_interp sources
#set(amip_interp_SOURCES amip_interp/amip_interp.F90)

# **********************************
# astronomy sources
#set(astronomy_SOURCES astronomy/astronomy.F90)
#list(APPEND libfms_SOURCES ${astronomy_SOURCES})

# **********************************
# axis_utils sources
#set(axis_utils_SOURCES axis_utils/axis_utils.F90)
#list(APPEND libfms_SOURCES ${axis_utils_SOURCES})

# **********************************
# block_control sources
#set(block_control_SOURCES block_control/block_control.F90)
#list(APPEND libfms_SOURCES ${block_control_SOURCES})

# **********************************
# column_diagnostics sources
#set(column_diagnostics_SOURCES column_diagnostics/column_diagnostics.F90)
#list(APPEND libfms_SOURCES ${column_diagnostics_SOURCES})

# **********************************
# constants sources
#set(constants_SOURCES constants/constants.F90)
#list(APPEND libfms_SOURCES ${constants_SOURCES})

# **********************************
# coupler sources
#set(coupler_SOURCES
#  coupler/atmos_ocean_fluxes.F90
#  coupler/coupler_types.F90
#  coupler/ensemble_manager.F90
#  )
#list(APPEND libfms_SOURCES ${coupler_SOURCES})

# **********************************
# data_override sources
#set(data_override_SOURCES data_override/data_override.F90)
#list(APPEND libfms_SOURCES ${data_override_SOURCES})

# **********************************
# diag_manager sources
#set(diag_manager_SOURCES
#  diag_manager/diag_axis.F90
#  diag_manager/diag_data.F90
#  diag_manager/diag_grid.F90
#  diag_manager/diag_manager.F90
#  diag_manager/diag_output.F90
#  diag_manager/diag_table.F90
#  diag_manager/diag_util.F90
#  diag_manager/diag_manifest.F90
#  )
#list(APPEND libfms_SOURCES ${diag_manager_SOURCES})

# **********************************
# drifters sources
#set(drifters_SOURCES
#  drifters/cloud_interpolator.F90
#  drifters/drifters.F90
#  drifters/drifters_comm.F90
#  drifters/drifters_core.F90
#  drifters/drifters_input.F90
#  drifters/drifters_io.F90
#  drifters/quicksort.F90
#  )
#list(APPEND libfms_SOURCES ${drifters_SOURCES})

# **********************************
# exchange sources
# set(exchange_SOURCES
#   exchange/stock_constants.F90
#   exchange/test_xgrid.F90
#   exchange/xgrid.F90
#   )
# list(APPEND libfms_SOURCES ${exchange_SOURCES})

# **********************************
# fft sources
# set(fft_SOURCES
#   fft/fft.F90
#   fft/fft99.F90
#   )
# list(APPEND libfms_SOURCES ${fft_SOURCES})

# **********************************
# field_manager sources
# include_directories("${PROJECT_SOURCE_DIR}/fms")
# set(field_manager_SOURCES
#   field_manager/field_manager.F90
#   field_manager/fm_util.F90
#   )
# list(APPEND libfms_SOURCES ${field_manager_SOURCES})

# **********************************
# fms sources
# set(fms_SOURCES
#   fms/fms.F90
#   fms/fms_io.F90
#   )
# list(APPEND libfms_SOURCES ${fms_SOURCES})

# **********************************
# horiz_interp sources
# set(horiz_interp_SOURCES
#   horiz_interp/horiz_interp.F90
#   horiz_interp/horiz_interp_bicubic.F90
#   horiz_interp/horiz_interp_bilinear.F90
#   horiz_interp/horiz_interp_conserve.F90
#   horiz_interp/horiz_interp_spherical.F90
#   horiz_interp/horiz_interp_type.F90
#   )
# list(APPEND libfms_SOURCES ${horiz_interp_SOURCES})

# **********************************
# interpolator sources
#set(interpolator_SOURCES interpolator/interpolator.F90)
#list(APPEND libfms_SOURCES ${interpolator_SOURCES})

# **********************************
# memutils sources
# set(memutils_SOURCES
#   memutils/memuse.c
#   memutils/memutils.F90
#   )
# list(APPEND libfms_SOURCES ${memutils_SOURCES})

# **********************************
# mosaic sources
# set(mosaic_SOURCES
#   mosaic/create_xgrid.c
#   mosaic/gradient.F90
#   mosaic/gradient_c2l.c
#   mosaic/grid.F90
#   mosaic/interp.c
#   mosaic/mosaic.F90
#   mosaic/mosaic_util.c
#   mosaic/read_mosaic.c
#   )
# list(APPEND libfms_SOURCES ${mosaic_SOURCES})

# **********************************
# mpp sources
# include_directories("${PROJECT_SOURCE_DIR}/mpp/include")
# set(mpp_SOURCES
#   mpp/affinity.c
#   mpp/nsclock.c
#   mpp/threadloc.c
#   mpp/mpp_domains.F90
#   mpp/mpp_data.F90
#   mpp/mpp_efp.F90
#   mpp/mpp.F90
#   mpp/mpp_io.F90
#   mpp/mpp_memutils.F90
#   mpp/mpp_parameter.F90
#   mpp/mpp_pset.F90
#   mpp/mpp_utilities.F90
#   )
# set(mpp_CPP "-Duse_libMPI -Duse_netCDF")

# list(APPEND libfms_SOURCES ${mpp_SOURCES})

# **********************************
# oda_tools sources
# set(oda_tools_SOURCES
#   oda_tools/oda_core.F90
#   oda_tools/oda_core_ecda.F90
#   oda_tools/oda_types.F90
#   oda_tools/write_ocean_data.F90
#   oda_tools/xbt_drop_rate_adjust.f90
#   )
# list(APPEND libfms_SOURCES ${oda_tools_SOURCES})

# **********************************
# platform sources
#set(platform_SOURCES platform/platform.F90)
#list(APPEND libfms_SOURCES ${platform_SOURCES})

# **********************************
# random_numbers sources
# set(random_numbers_SOURCES
#   random_numbers/MersenneTwister.F90
#   random_numbers/random_numbers.F90
#   )
# list(APPEND libfms_SOURCES ${random_numbers_SOURCES})

# **********************************
# sat_vapor_pres sources
# set(sat_vapor_pres_SOURCES
#   sat_vapor_pres/sat_vapor_pres.F90
#   sat_vapor_pres/sat_vapor_pres_k.F90
#   )
# list(APPEND libfms_SOURCES ${sat_vapor_pres_SOURCES})

# **********************************
# time_interp sources
# set(time_interp_SOURCES
#   time_interp/time_interp.F90
#   time_interp/time_interp_external.F90
#   )
# list(APPEND libfms_SOURCES ${time_interp_SOURCES})

# **********************************
# time_manager sources
# set(time_manager_SOURCES
#   time_manager/get_cal_time.F90
#   time_manager/time_manager.F90
#   )
# list(APPEND libfms_SOURCES ${time_manager_SOURCES})

# **********************************
# topography sources
# set(topography_SOURCES
#   topography/gaussian_topog.F90
#   topography/topography.F90
#   )
# list(APPEND libfms_SOURCES ${topography_SOURCES})

# **********************************
# tracer_manager sources
#set(tracer_manager_SOURCES tracer_manager/tracer_manager.F90)
#list(APPEND libfms_SOURCES ${tracer_manager_SOURCES})

# **********************************
# tridiagonal sources
#set(tridiagonal_SOURCES tridiagonal/tridiagonal.F90)
#list(APPEND libfms_SOURCES ${tridiagonal_SOURCES})

# **********************************************************************
# Add library to package

# **********************************************************************
# Set Compiler/Linker Flags

#set(CMAKE_Fortran_FLAGS "${GNU_F90_FLAGS_BASE}")

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
