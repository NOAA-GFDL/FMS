project(time_manager Fortran)

set(${PROJECT_NAME}_SOURCES
  get_cal_time.F90
  time_manager.F90
  )

include_directories(
  ${constants_INCLUDE_PATH}
  ${mpp_INCLUDE_PATH}
  ${fms_INCLUDE_PATH}
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms)

# Expose public includes to other subprojects through a cache variable.
set(${PROJECT_NAME}_INCLUDE_PATH ${PROJECT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
