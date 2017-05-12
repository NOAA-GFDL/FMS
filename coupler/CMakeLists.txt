project(coupler Fortran)

set(${PROJECT_NAME}_SOURCES
  atmos_ocean_fluxes.F90
  coupler_types.F90
  ensemble_manager.F90
  )

include_directories(
  ${mpp_INCLUDE_PATH}
  ${fms_INCLUDE_PATH}
  ${field_manager_INCLUDE_PATH}
  ${time_manager_INCLUDE_PATH}
  ${diag_manager_INCLUDE_PATH}
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} mpp fms field_manager time_manager diag_manager)

# Expose public includes to other subprojects through a cache variable.
set(${PROJECT_NAME}_INCLUDE_PATH ${PROJECT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
