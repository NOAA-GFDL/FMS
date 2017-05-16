project(topography Fortran)

set(${PROJECT_NAME}_SOURCES
  gaussian_topog.F90
  topography.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms_mod horiz_interp)

