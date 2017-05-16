project(astronomy Fortran)

set(${PROJECT_NAME}_SOURCES
  astronomy.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms_mod time_manager)

