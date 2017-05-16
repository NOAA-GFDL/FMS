project(oda_tools Fortran)

set(${PROJECT_NAME}_SOURCES
  oda_core_ecda.F90
  oda_core.F90
  oda_types.F90
  write_ocean_data.F90
  xbt_drop_rate_adjust.f90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms_mod axis_utils field_manager time_manager horiz_interp)

if (NETCDF_Fortran_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC "use_netCDF")
  include_directories(${NETCDF_Fortran_INCLUDE_PATH})
endif (NETCDF_Fortran_FOUND)

