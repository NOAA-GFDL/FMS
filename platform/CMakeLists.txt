project(platform Fortran)

set(${PROJECT_NAME}_SOURCES platform.F90)

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})

# Expose public includes to other subprojects through a cache variable.
set(${PROJECT_NAME}_INCLUDE_PATH ${PROJECT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
