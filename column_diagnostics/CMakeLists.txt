project(column_diagnostics Fortran)

set(${PROJECT_NAME}_SOURCES
  column_diagnostics.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms_mod time_manager)

