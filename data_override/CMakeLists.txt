project(data_override Fortran)

set(${PROJECT_NAME}_SOURCES
  data_override.F90
  )

include_directories(
  ${constants_INCLUDE_PATH}
  ${platform_INCLUDE_PATH}
  ${mpp_INCLUDE_PATH}
  ${fms_INCLUDE_PATH}
  ${axis_utils_INCLUDE_PATH}
  ${time_manager_INCLUDE_PATH}
  ${diag_manager_INCLUDE_PATH}
  ${horiz_interp_INCLUDE_PATH}
  ${time_interp_INCLUDE_PATH}
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants platform mpp fms axis_utils time_manager diag_manager horiz_interp time_interp)

# Expose public includes to other subprojects through a cache variable.
set(${PROJECT_NAME}_INCLUDE_PATH ${PROJECT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
