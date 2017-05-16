project(sat_vapor_pres Fortran)

set(${PROJECT_NAME}_SOURCES
  sat_vapor_pres.F90
  sat_vapor_pres_k.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms_mod)

