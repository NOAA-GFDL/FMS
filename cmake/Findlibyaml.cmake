# Try to find libyaml includes and library.
# The following variables are set:
#    LIBYAML_INCLUDE_DIR
#    LIBYAML_LIBRARIES

FIND_PATH(LIBYAML_INCLUDE_DIR  NAMES yaml.h  PATHS ${LIBYAML_ROOT}/include $ENV{LIBYAML_ROOT}/include )
FIND_LIBRARY(LIBYAML_LIBRARIES NAMES yaml    PATHS ${LIBYAML_ROOT}/lib     $ENV{LIBYAML_ROOT}/lib )
if(NOT LIBYAML_INCLUDE_DIR OR NOT LIBYAML_LIBRARIES)
  message(SEND_ERROR "libyaml library/include file not found, set LIBYAML_ROOT")
endif()
MARK_AS_ADVANCED(LIBYAML_INCLUDE_DIR LIBYAML_LIBRARIES)
