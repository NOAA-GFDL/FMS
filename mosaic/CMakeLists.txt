project(mosaic C Fortran)

set(${PROJECT_NAME}_SOURCES
  constant.h
  create_xgrid.c
  create_xgrid.h
  gradient_c2l.c
  gradient_c2l.h
  gradient.F90
  grid.F90
  interp.c
  interp.h
  mosaic.F90
  mosaic_util.c
  mosaic_util.h
  read_mosaic.c
  read_mosaic.h
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms_mod)

if (NETCDF_Fortran_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC "use_netCDF")
  include_directories(${NETCDF_Fortran_INCLUDE_PATH})
endif (NETCDF_Fortran_FOUND)

if (MPI_Fortran_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC "use_MPI")
  include_directories(${MPI_Fortran_INCLUDE_PATH})
endif (MPI_Fortran_FOUND)

