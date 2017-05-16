project(tracer_manager Fortran)

set(${PROJECT_NAME}_SOURCES
  tracer_manager.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} mpp fms_mod field_manager)

