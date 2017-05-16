project(diag_manager Fortran)

set(${PROJECT_NAME}_SOURCES
  diag_axis.F90
  diag_data.F90
  diag_grid.F90
  diag_manager.F90
  diag_manifest.F90
  diag_output.F90
  diag_table.F90
  diag_util.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms_mod time_manager)

if (NETCDF_Fortran_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC "use_netCDF")
  include_directories(${NETCDF_Fortran_INCLUDE_PATH})
endif (NETCDF_Fortran_FOUND)
