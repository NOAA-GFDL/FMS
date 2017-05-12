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

include_directories(
  ${constants_INCLUDE_PATH}
  ${mpp_INCLUDE_PATH}
  ${fms_INCLUDE_PATH}
  )

add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})
add_dependencies(${PROJECT_NAME} constants mpp fms)

if (NETCDF_FOUND OR NETCDFF_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC "use_netCDF")
  include_directories(${NETCDF_INCLUDE_PATH} ${NETCDFF_INCLUDE_PATH})
endif (NETCDF_FOUND OR NETCDFF_FOUND)

if (MPI_C_FOUND OR MPI_Fortran_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC "use_MPI")
  include_directories(${MPI_C_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
endif (MPI_C_FOUND OR MPI_Fortran_FOUND)

# Expose public includes to other subprojects through a cache variable.
set(${PROJECT_NAME}_INCLUDE_PATH ${PROJECT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
