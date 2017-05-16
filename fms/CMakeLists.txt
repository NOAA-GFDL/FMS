project(fms_mod C Fortran)

set(${PROJECT_NAME}_SOURCES
  fms.F90
  fms_io.F90
  )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants platform memutils mpp)
