#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "FMS::fms_r4" for configuration "Release"
set_property(TARGET FMS::fms_r4 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(FMS::fms_r4 PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libfms_r4.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS FMS::fms_r4 )
list(APPEND _IMPORT_CHECK_FILES_FOR_FMS::fms_r4 "${_IMPORT_PREFIX}/lib/libfms_r4.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
