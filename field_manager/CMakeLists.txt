project(field_manager Fortran)

set(${PROJECT_NAME}_SOURCES
  field_manager.F90
  fm_util.F90
  parse.inc
  )

include_directories(${mpp_INCLUDE_PATH}
  ${fms_INCLUDE_PATH}
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} mpp fms)

# Expose public includes to other subprojects through a cache variable.
set(${PROJECT_NAME}_INCLUDE_PATH ${PROJECT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
