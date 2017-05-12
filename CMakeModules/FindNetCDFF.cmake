find_program(NF_CONFIG
  NAMES nf-config
  DOC "Executable for finding NetCDF Fortran configuration options")

if (NETCDFF_INCLUDE_PATH AND NETCDFF_LIBRARY)
  # Do nothing: we already have NETCDF_INCLUDE_PATH and NETCDF_LIBRARY in
  # the cache, and we don't want to override those settings.
elseif (NF_CONFIG)
  exec_program(${NF_CONFIG}
    ARGS --includedir
    OUTPUT_VARIABLE NETCDFF_INCLUDE_PATH
    RETURN_VALUE NETCDFF_INCLUDE_RETURN)

  exec_program(${NF_CONFIG}
    ARGS --flibs
    OUTPUT_VARIABLE NETCDFF_LIBRARY
    RETURN_VALUE NETCDFF_LIBRARY_RETURN)
endif (NETCDFF_INCLUDE_PATH AND NETCDFF_LIBRARY)

if (NETCDFF_INCLUDE_PATH AND NETCDFF_LIBRARY)
  set(NETCDFF_FOUND TRUE)
else (NETCDFF_INCLUDE_PATH AND NETCDFF_LIBRARY)
  set(NETCDFF_FOUND FALSE)
endif (NETCDFF_INCLUDE_PATH AND NETCDFF_LIBRARY)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments 
find_package_handle_standard_args(NetCDFF DEFAULT_MSG NETCDFF_LIBRARY NETCDFF_INCLUDE_PATH)

mark_as_advanced(NETCDFF_LIBRARY NETCDFF_INCLUDE_PATH)
