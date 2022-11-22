
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was FMSConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

#fms-config.cmake
#
# Valid Find COMPONENTS:
#  * R4 - 32-bit library
#  * R8 - 64-bit library
#
# Imported interface targets provided:
#  * FMS::fms_r4 - real32 library target
#  * FMS::fms_r8 - real64 library target

# Include targets file.  This will create IMPORTED target fms
include("${CMAKE_CURRENT_LIST_DIR}/fms-targets.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/fms-config-version.cmake")
include(CMakeFindDependencyMacro)

find_dependency(MPI COMPONENTS Fortran)

# ON/OFF implies FMS was compiled with/without OpenMP
if(OFF)
  find_dependency(OpenMP COMPONENTS Fortran)
endif()

find_dependency(NetCDF COMPONENTS C Fortran)

set(FMSVersion "${PACKAGE_VERSION}")
set_and_check(FMS_INSTALL_PREFIX "${PACKAGE_PREFIX_DIR}")

list(APPEND _possible_fms_components "R4" "R8")
foreach(_comp ${_possible_fms_components})
  set(_name_${_comp} ${_comp})
endforeach()

#R4 is the default component
if(NOT FMS_FIND_COMPONENTS)
  list(APPEND _search_fms_components "R4")
endif()
foreach(_comp ${FMS_FIND_COMPONENTS})
  list(APPEND _search_fms_components ${_comp})
  if(NOT _name_${_comp})
    message(SEND_ERROR "FMS: COMPONENT ${_comp} is not a valid component. Valid components: ${_possible_fms_components}")
  endif()
endforeach()
list(REMOVE_DUPLICATES _search_fms_components)

# ON/OFF implies FMS was compiled with/without 32BIT option
set(FMS_R4_FOUND 0)
if(ON)
  set(FMS_R4_FOUND 1)
  get_property(FMS_R4_LIBRARIES TARGET FMS::fms_r4 PROPERTY LOCATION)
  get_target_property(FMS_BUILD_TYPES FMS::fms_r4 IMPORTED_CONFIGURATIONS)
  list(APPEND FMS_POSSIBLE_COMPONENTS R4)
endif()

# ON/OFF implies FMS was compiled with/without 64BIT option
set(FMS_R8_FOUND 0)
if(OFF)
  set(FMS_R8_FOUND 1)
  get_property(FMS_R8_LIBRARIES TARGET FMS::fms_r8 PROPERTY LOCATION)
  get_target_property(FMS_BUILD_TYPES FMS::fms_r8 IMPORTED_CONFIGURATIONS)
endif()

check_required_components("FMS")

message(STATUS "Found FMS: \"${FMS_INSTALL_PREFIX}\" (Version: \"${FMSVersion}\")")
message( STATUS "FMS targets:" )
foreach(component ${_search_fms_components})
  string( TOLOWER "${component}" _component )
  if(FMS_${component}_FOUND)
    message(STATUS "  - FMS::fms_${_component} [Lib: ${FMS_${component}_LIBRARIES}]")
  else()
    message(STATUS "  - FMS::fms_${_component} [COMPONENT ${component} NOT FOUND]")
  endif()
endforeach()
