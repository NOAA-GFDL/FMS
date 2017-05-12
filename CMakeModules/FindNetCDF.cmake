find_program(NC_CONFIG
  NAMES nc-config
  DOC "Executable for finding NetCDF configuration options")

if (NETCDF_INCLUDE_PATH AND NETCDF_LIBRARY)
  # Do nothing: we already have NETCDF_INCLUDE_PATH and NETCDF_LIBRARY in
  # the cache, and we don't want to override those settings.
elseif (NC_CONFIG)
  exec_program(${NC_CONFIG}
    ARGS --includedir
    OUTPUT_VARIABLE NETCDF_INCLUDE_PATH
    RETURN_VALUE NETCDF_INCLUDE_RETURN)

  exec_program(${NC_CONFIG}
    ARGS --flibs
    OUTPUT_VARIABLE NETCDF_LIBRARY
    RETURN_VALUE NETCDF_LIBRARY_RETURN)
endif (NETCDF_INCLUDE_PATH AND NETCDF_LIBRARY)

if (NETCDF_INCLUDE_PATH AND NETCDF_LIBRARY)
  set(NETCDF_FOUND TRUE)
else (NETCDF_INCLUDE_PATH AND NETCDF_LIBRARY)
  set(NETCDF_FOUND FALSE)
endif (NETCDF_INCLUDE_PATH AND NETCDF_LIBRARY)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments 
find_package_handle_standard_args(NetCDF DEFAULT_MSG NETCDF_LIBRARY NETCDF_INCLUDE_PATH)

mark_as_advanced(NETCDF_LIBRARY NETCDF_INCLUDE_PATH)