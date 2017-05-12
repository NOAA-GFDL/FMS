project(drifters Fortran)

set(${PROJECT_NAME}_SOURCES
  cloud_interpolator.F90
  drifters.F90
  drifters_comm.F90
  drifters_core.F90
  drifters_input.F90
  drifters_io.F90
  quicksort.F90
  )

include_directories(${mpp_INCLUDE_PATH})

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} mpp)

if (NETCDF_FOUND OR NETCDFF_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC "use_netCDF")
  include_directories(${NETCDF_INCLUDE_PATH} ${NETCDFF_INCLUDE_PATH})
endif (NETCDF_FOUND OR NETCDFF_FOUND)

# Expose public includes to other subprojects through a cache variable.
set(${PROJECT_NAME}_INCLUDE_PATH ${PROJECT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
