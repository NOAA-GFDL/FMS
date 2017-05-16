project(time_manager Fortran)

set(${PROJECT_NAME}_SOURCES
  get_cal_time.F90
  time_manager.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms_mod)

