project(memutils C Fortran)

set(${PROJECT_NAME}_SOURCES memuse.c memutils.F90)

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} mpp)

