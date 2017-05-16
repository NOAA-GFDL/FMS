project(data_override Fortran)

set(${PROJECT_NAME}_SOURCES
  data_override.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants platform mpp fms_mod axis_utils time_manager diag_manager horiz_interp time_interp)

