project(time_interp Fortran)

set(${PROJECT_NAME}_SOURCES
  time_interp_external.F90
  time_interp.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants platform mpp fms_mod axis_utils time_manager horiz_interp)

