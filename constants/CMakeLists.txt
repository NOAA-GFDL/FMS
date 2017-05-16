project(constants Fortran)

set(${PROJECT_NAME}_SOURCES constants.F90)

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})

