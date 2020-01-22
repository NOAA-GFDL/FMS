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

cmake_minimum_required(VERSION 3.7..3.15 FATAL_ERROR)

# Define the CMake project
project(libFMS
  VERSION 2020.02
  DESCRIPTION  "GFDL FMS Library"
  HOMEPAGE_URL https://www.gfdl.noaa.gov/fms
  LANGUAGES C Fortran)

add_definitions(-Duse_libMPI)
add_definitions(-Duse_netCDF)

# Set the build options
# These can be changed via `cmake -D<var>=[ON|OFF]`, e.g. cmake -DOPENMP=OFF to
# disable OpenMP
option(USE_OPENMP "Determine if libFMS should be built with OpenMP support" ON)

# Add directory that contains custom CMake Modules for building libFMS
list(INSERT CMAKE_MODULE_PATH 0 ${libFMS_SOURCE_DIR}/CMakeModules)

# Check for the MPI libaray, and compiler wrappers
find_package(MPI REQUIRED)

# Check for the OpenMP library and set the required compile flags
if (OPENMP)
  find_package(OpenMP REQUIRED)
  if (OPENMP_FOUND)
    list(APPEND CMAKE_C_FLAGS ${OpenMP_C_FLAGS})
    list(APPEND CMAKE_Fortran_FLAGS ${OpenMP_Fortran_FLAGS})
  endif (OPENMP_FOUND)
endif(OPENMP)

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${MPI_C_INCLUDE_PATH}
  ${MPI_Fortran_INCLUDE_PATH}
)

message("************************************************************")
message("Building FMS with")
message("Fortran Compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message("Fortran Flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}")
message("C Compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message("C Flags: ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}")
message("Build Target: ${BUILD_TARGET}")
message("OpenMP Enabled: ${OPENMP}")
message("Using NetCDF C: ${NETCDF_C_VERSION}")
message("Using NetCDF Fortran: ${NETCDF_Fortran_VERSION}")
message("MPI: ${MPI_C_FOUND} ${MPI_Fortran_FOUND}")
message("MPI Compiler: ${MPI_C_COMPILER} ${MPI_Fortran_COMPILER}")
message("MPI Compiler Flags: ${MPI_C_COMPILE_FLAGS} | ${MPI_F_COMPILE_FLAGS}")
message("MPI INCLUDES: ${MPI_C_INCLUDE_PATH} | ${MPI_Fortran_INCLUDE_PATH}")
message("************************************************************")

add_library(FMS
  affinity/affinity.c
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
  diag_manager/diag_manifest.F90
  diag_manager/diag_output.F90
  diag_manager/diag_table.F90
  diag_manager/diag_util.F90
  exchange/stock_constants.F90
  exchange/xgrid.F90
  field_manager/field_manager.F90
  field_manager/fm_util.F90
  fms/fms_io.F90
  fms/fms.F90
  fms2_io/blackboxio.f90
  fms2_io/fms_io_utils.F90
  fms2_io/fms_netcdf_domain_io.f90
  fms2_io/fms_netcdf_unstructured_domain_io.f90
  fms2_io/fms2_io.f90
  fms2_io/legacy.f90
  fms2_io/netcdf_io.f90
  horiz_interp/horiz_interp_bicubic.F90
  horiz_interp/horiz_interp_bilinear.F90
  horiz_interp/horiz_interp_conserve.F90
  horiz_interp/horiz_interp_spherical.F90
  horiz_interp/horiz_interp_type.f90
  horiz_interp/horiz_interp.F90
  interpolator/interpolator.F90
  memutils/memuse.c
  memutils/memutils.F90
  monin_obukhov/monin_obukhov_kernel.F90
  monin_obukhov/monin_obukhov.f90
  mosaic/create_xgrid.c
  mosaic/gradient_c2l.c
  mosaic/gradient.F90
  mosaic/grid.F90
  mosaic/interp.c
  mosaic/mosaic_util.c
  mosaic/mosaic.F90
  mosaic/mosaic2.F90
  mosaic/read_mosaic.c
  mpp/mpp.F90
  mpp/mpp_data.F90
  mpp/mpp_domains.F90
  mpp/mpp_efp.F90
  mpp/mpp_io.F90
  mpp/mpp_memutils.F90
  mpp/mpp_parameter.F90
  mpp/mpp_pset.F90
  mpp/mpp_utilities.F90
  mpp/nsclock.c
  mpp/threadloc.c
  platform/platform.F90
  random_numbers/MersenneTwister.f90
  random_numbers/random_numbers.f90
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
  tridiagonal/tridiagonal.f90
)
