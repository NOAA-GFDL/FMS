project(platform Fortran)

set(${PROJECT_NAME}_SOURCES platform.F90)

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})

