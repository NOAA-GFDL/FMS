project(exchange Fortran)

set(${PROJECT_NAME}_SOURCES
  stock_constants.F90
  test_xgrid.F90
  xgrid.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms_mod mosaic time_manager diag_manager coupler)
