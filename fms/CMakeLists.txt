project(fms C Fortran)

set(${PROJECT_NAME}_SOURCES
  fms.F90
  fms_io.F90
  )

include_directories(${mpp_INCLUDE_PATH}
  ${constants_INCLUDE_PATH}
  ${platform_INCLUDE_PATH}
  ${memutils_INCLUDE_PATH}
  )
add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants platform memutils mpp)

# Expose public includes to other subprojects through a cache variable.
set(${PROJECT_NAME}_INCLUDE_PATH ${PROJECT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
