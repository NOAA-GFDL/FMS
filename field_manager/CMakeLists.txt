project(field_manager Fortran)

set(${PROJECT_NAME}_SOURCES
  field_manager.F90
  fm_util.F90
  parse.inc
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms_mod)

