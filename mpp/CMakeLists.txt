project(mpp C Fortran)

set(${PROJECT_NAME}_SOURCES
  affinity.c
  nsclock.c
  threadloc.c
  mpp_domains.F90
  mpp_data.F90
  mpp_efp.F90
  mpp.F90
  mpp_io.F90
  mpp_memutils.F90
  mpp_parameter.F90
  mpp_pset.F90
  mpp_utilities.F90
  )
  
add_library(${PROJECT_NAME} OBJECT ${${PROJECT_NAME}_SOURCES})

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
