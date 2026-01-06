# Try to find libyaml includes and library.
# The following variables are set:
#   libyaml_FOUND
#   LIBYAML_INCLUDE_DIRS
#   LIBYAML_LIBRARIES
#
# Imported targets:
#   libyaml::yaml
#
# Hints:
#   LIBYAML_ROOT or ENV{LIBYAML_ROOT}

include(FindPackageHandleStandardArgs)

# Allow user to hint the install prefix
set(_LIBYAML_HINTS
  ${LIBYAML_ROOT}
  $ENV{LIBYAML_ROOT}
)

# Find the header
find_path(LIBYAML_INCLUDE_DIR
  NAMES yaml.h
  HINTS ${_LIBYAML_HINTS}
  PATH_SUFFIXES include
)

# Find the library
find_library(LIBYAML_LIBRARIES
  NAMES yaml
  HINTS ${_LIBYAML_HINTS}
  PATH_SUFFIXES lib lib64
)

# Standard handling of REQUIRED / QUIET and *_FOUND
find_package_handle_standard_args(libyaml
  REQUIRED_VARS
    LIBYAML_INCLUDE_DIR
    LIBYAML_LIBRARIES
)

if(libyaml_FOUND)
  if(NOT TARGET libyaml::yaml)
    add_library(libyaml::yaml UNKNOWN IMPORTED)
    set_target_properties(libyaml::yaml PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES ${LIBYAML_INCLUDE_DIR}
      IMPORTED_LOCATION             ${LIBYAML_LIBRARIES}
    )
  endif()
endif()

mark_as_advanced(
  LIBYAML_INCLUDE_DIR
  LIBYAML_LIBRARIES
)
