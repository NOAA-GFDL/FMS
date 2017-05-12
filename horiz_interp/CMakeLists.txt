project(horiz_interp Fortran)

set(${PROJECT_NAME}_SOURCES
  horiz_interp_bicubic.F90
  horiz_interp_bilinear.F90
  horiz_interp_conserve.F90
  horiz_interp.F90
  horiz_interp_spherical.F90
  horiz_interp_type.F90
  )

include_directories(
  ${constants_INCLUDE_PATH}
  ${mpp_INCLUDE_PATH}
  ${fms_INCLUDE_PATH}
  ${axis_utils_INCLUDE_PATH}
  ${mosaic_INCLUDE_PATH}
  ${time_manager_INCLUDE_PATH}
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms axis_utils mosaic time_manager)

# Expose public includes to other subprojects through a cache variable.
set(${PROJECT_NAME}_INCLUDE_PATH ${PROJECT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
