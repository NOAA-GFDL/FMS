project(drifters Fortran)

set(${PROJECT_NAME}_SOURCES
  cloud_interpolator.F90
  drifters.F90
  drifters_comm.F90
  drifters_core.F90
  drifters_input.F90
  drifters_io.F90
  quicksort.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} mpp)

if (NETCDF_Fortran_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC "use_netCDF")
  include_directories(${NETCDF_Fortran_INCLUDE_PATH})
endif (NETCDF_Fortran_FOUND)
