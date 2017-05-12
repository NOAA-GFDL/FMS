project(oda_tools Fortran)

set(${PROJECT_NAME}_SOURCES
  oda_core_ecda.F90
  oda_core.F90
  oda_types.F90
  write_ocean_data.F90
  xbt_drop_rate_adjust.f90
  )

include_directories(
  ${constants_INCLUDE_PATH}
  ${mpp_INCLUDE_PATH}
  ${fms_INCLUDE_PATH}
  ${axis_utils_INCLUDE_PATH}
  ${field_manager_INCLUDE_PATH}
  ${time_manager_INCLUDE_PATH}
  ${horiz_interp_INCLUDE_PATH}
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms axis_utils field_manager time_manager horiz_interp)

if (NETCDF_FOUND OR NETCDFF_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC "use_netCDF")
  include_directories(${NETCDF_INCLUDE_PATH} ${NETCDFF_INCLUDE_PATH})
endif (NETCDF_FOUND OR NETCDFF_FOUND)

# Expose public includes to other subprojects through a cache variable.
set(${PROJECT_NAME}_INCLUDE_PATH ${PROJECT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
