project(exchange Fortran)

set(${PROJECT_NAME}_SOURCES
  stock_constants.F90
  test_xgrid.F90
  xgrid.F90
  )

include_directories(
  ${constants_INCLUDE_PATH}
  ${mpp_INCLUDE_PATH}
  ${fms_INCLUDE_PATH}
  ${mosaic_INCLUDE_PATH}
  ${time_manager_INCLUDE_PATH}
  ${diag_manager_INCLUDE_PATH}
  ${coupler_INCLUDE_PATH}
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms mosaic time_manager diag_manager coupler)

# Expose public includes to other subprojects through a cache variable.
set(${PROJECT_NAME}_INCLUDE_PATH ${PROJECT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
