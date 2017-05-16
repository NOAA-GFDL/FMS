project(coupler Fortran)

set(${PROJECT_NAME}_SOURCES
  atmos_ocean_fluxes.F90
  coupler_types.F90
  ensemble_manager.F90
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} mpp fms_mod field_manager time_manager diag_manager)

