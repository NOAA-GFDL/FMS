project(diag_manager Fortran)

set(${PROJECT_NAME}_SOURCES
  diag_axis.F90
  diag_data.F90
  diag_grid.F90
  diag_manager.F90
  diag_manifest.F90
  diag_output.F90
  diag_table.F90
  diag_util.F90
  )

include_directories(
  ${constants_INCLUDE_PATH}
  ${mpp_INCLUDE_PATH}
  ${fms_INCLUDE_PATH}
  ${time_manager_INCLUDE_PATH}
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms time_manager)

if (NETCDF_FOUND OR NETCDFF_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC "use_netCDF")
  include_directories(${NETCDF_INCLUDE_PATH} ${NETCDFF_INCLUDE_PATH})
endif (NETCDF_FOUND OR NETCDFF_FOUND)

# Expose public includes to other subprojects through a cache variable.
set(${PROJECT_NAME}_INCLUDE_PATH ${PROJECT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
