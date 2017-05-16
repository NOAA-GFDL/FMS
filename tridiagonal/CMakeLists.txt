project(tridiagonal Fortran)

set(${PROJECT_NAME}_SOURCES
  tridiagonal.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})

