project(block_control Fortran)

set(${PROJECT_NAME}_SOURCES block_control.F90)

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} mpp)

