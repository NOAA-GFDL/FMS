project(random_numbers Fortran)

set(${PROJECT_NAME}_SOURCES
  MersenneTwister.F90
  random_numbers.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} time_manager)

