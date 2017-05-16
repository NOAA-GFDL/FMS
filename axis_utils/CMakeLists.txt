project(axis_utils Fortran)

set(${PROJECT_NAME}_SOURCES
  axis_utils.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms_mod)

if (NETCDF_Fortran_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC "use_netCDF")
  include_directories(${NETCDF_Fortran_INCLUDE_PATH})
endif (NETCDF_Fortran_FOUND)
