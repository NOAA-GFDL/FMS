project(random_numbers Fortran)

set(${PROJECT_NAME}_SOURCES
  MersenneTwister.F90
  random_numbers.F90
  )

include_directories(
  ${time_manager_INCLUDE_PATH}
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} time_manager)

# Expose public includes to other subprojects through a cache variable.
set(${PROJECT_NAME}_INCLUDE_PATH ${PROJECT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
