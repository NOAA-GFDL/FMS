project(amip_interp Fortran)

set(${PROJECT_NAME}_SOURCES
  amip_interp.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants platform mpp fms_mod time_manager horiz_interp time_interp)

