# Distributed with the Flexible Modeling System (FMS) source code
# Author: Seth Underwood

#.rst
# FindNetCDF
# ----------
#
# Find the NetCDF libraries
#
# The Network Common Data Form (netCDF) is a set of software libraries
# and self-describing, machine-independent data formats that support
# the creation, access, and sharing of array-oriented scientific
# data. NetCDF is maintained by Unidata program at the University
# Corporation for Atmospheric Research (UCAR), and is obtainable at
# https://www.unidata.ucar.edu/software/netcdf/
#
# Variables
# ^^^^^^^^^
#
# This module will set the following variables per language in your
# project, where <lang> is one of C or Fortran.
#
# ::
#
#    NETCDF_<lang>_FOUND            TRUE if FindNetCDF found netCDF for <lang>
#    NETCDF_<lang>_COMPILER         Compiler used to build netCDF for <lang>
#    NETCDF_<lang>_COMPILE_FLAGS    Compilation flags for netCDF
#    NETCDF_<lang>_INCLUDE_PATH     Include path(s) for netCDF headers for <lang>
#    NETCDF_<lang>_LIBRARIES        Link flags and all libraries to link netCDF programs against
#
# Usage
# ^^^^^
#
# To use this module, call FindNetCDF from a CMakeList.txt file, or
# run ``find_package(NetCDF)``.  This module requires the programs
# ``nc-confg`` and ``nf-config`` to automatically set the above
# variables for <lang> C and Fortran respecitvely.  If the programs
# cannot be found by cmake, you will need to set the variables
# manually:
#
# ::
#
#    NETCDF_<lang>_INCLUDE_PATH
#    NETCDF_<lang>_LIBRARIES
#
# If that still does not work, the user will need to also set
# NETCDF_<lang>_COMPILE_FLAGS.

# include this to handle the QUIETLY and REQUIRED arguments
include(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

find_program(NCDF_CONFIG_C
  NAMES nc-config
  HINTS ${NETCDF}/bin
  DOC "Executable for finding NetCDF configuration options")

find_program(NCDF_CONFIG_Fortran
  NAMES nf-config
  HINTS ${NETCDF}/bin
  DOC "Executable for finding NetCDF Fortran configuration options")

# If we could not find the netCDF configuration program, then the user
# _MUST_ set both NETCDF_<lang>_INCLUDE_PATH and NETCDF_<lang>_LIBRARIES
foreach (lang C Fortran)
  # Default REQUIRED flags set:
  if (NOT NCDF_CONFIG_${lang})
    if (NOT (NETCDF_${lang}_INCLUDE_PATH OR NETCDF_${lang}_LIBRARIES))
      if (NOT NETCDF_${lang}_COMPILER)
	set(NETCDF_${lang}_COMPILER CMAKE_${lang}_COMPILER)
      endif (NOT NETCDF_${lang}_COMPILER)
      if (NOT NETCDF_${lang}_COMPILE_FLAGS)
	set(NETCDF_${lang}_COMPILE_FLAGS "-I${NETCDF_${lang}_INCLUDE_PATH}")
      endif (NOT NETCDF_${lang}_COMPILE_FLAGS)
    endif (NOT (NETCDF_${lang}_INCLUDE_PATH OR NETCDF_${lang}_LIBRARIES))
  else (NOT NC_CONFIG_${lang})
    message("Doing language ${lang} loop")
    if (lang STREQUAL C)
      set(_libs_arg --libs)
      set(_comp_arg --cc)
      set(_flag_arg --cflags)
    elseif (lang STREQUAL Fortran)
      set(_libs_arg --flibs)
      set(_comp_arg --fc)
      set(_flag_arg --fflags)
    endif (lang STREQUAL C)
    if (NOT NETCDF_${lang}_INCLUDE_PATH)
      exec_program(${NCDF_CONFIG_${lang}}
	ARGS --includedir
	OUTPUT_VARIABLE NETCDF_INCLUDE_PATH
	RETURN_VALUE NETCDF_INCLUDE_RETURN
	)
    endif (NOT NETCDF_${lang}_INCLUDE_PATH)
    if (NOT NETCDF_${lang}_LIBRARIES)
      exec_program(${NCDF_CONFIG_${lang}}
	ARGS ${_libs_arg}
	OUTPUT_VARIABLE NETCDF_LIBRARIES
	RETURN_VALUE NETCDF_LIBRARIES_RETURN
	)
    endif (NOT NETCDF_${lang}_LIBRARIES)
    if (NOT NETCDF_${lang}_COMPILER)
      exec_program(${NCDF_CONFIG_${lang}}
	ARGS ${_comp_arg}
	OUTPUT_VARIABLE NETCDF_COMPILER
	RETURN_VALUE NETCDF_COMPILER_RETURN
	)
    endif (NOT NETCDF_${lang}_COMPILER)
    if (NOT NETCDF_${lang}_COMPILE_FLAGS)
      exec_program(${NCDF_CONFIG_${lang}}
	ARGS ${_flag_arg}
	OUTPUT_VARIABLE NETCDF_COMPILE_FLAGS
	RETURN_VALUE NETCDF_COMPILE_FLAGS_RETURN
	)
    endif (NOT NETCDF_${lang}_COMPILE_FLAGS)
    unset(_libs_arg)
    unset(_comp_arg)
    unset(_flag_arg)
    
    # Check the return value.  If either is zero, then C netCDF library
    # not found
    if (NOT (NETCDF_INCLUDE_RETURN EQUAL 0) OR
	NOT (NETCDF_LIBRARIES_RETURN EQUAL 0) OR
	NOT (NETCDF_COMPILER_RETURN EQUAL 0) OR
	NOT (NETCDF_COMPILE_FLAGS_RETURN EQUAL 0))
      set(NETCDF_${lang}_FOUND)
    else (NOT (NETCDF_INCLUDE_RETURN EQUAL 0) OR
	NOT (NETCDF_LIBRARIES_RETURN EQUAL 0) OR
	NOT (NETCDF_COMPILER_RETURN EQUAL 0) OR
	NOT (NETCDF_COMPILE_FLAGS_RETURN EQUAL 0))
      # Test the netCDF library
      set(scratch_directory ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY})
      if (lang STREQUAL C)
	set(test_file ${scratch_directory}/cmake_netcdf_test.c)
	file(WRITE ${test_file}
	  "#include <netcdf.h>\n"
	  "int main(int argc, char **argv) {\n"
	  "  char *fname;\n"
	  "  int *ncid;\n"
	  "  int ret;\n"
	  "  ret = nc_open(fname, NC_NOWRITE, ncid);\n"
	  "}\n")
      elseif (lang STREQUAL Fortran)
	set(test_file ${scratch_directory}/cmake_netcdf_test.F90)
	file(WRITE ${test_file}
	  "integer function f77_test()\n"
	  "#include <netcdf.inc>\n"
	  "character(len=32) :: fname\n"
	  "integer :: ncid\n"
	  "f77_test = nf_open(fname, NF_NOWRITE, ncid)\n"
	  "end function f77_test\n"
	  "integer function f90_test()\n"
	  "use netcdf\n"
	  "character(len=32) :: fname\n"
	  "integer :: ncid\n"
	  "f90_test = nf90_open(fname, NF90_NOWRITE, ncid)\n"
	  "end function f90_test\n"
	  "program test\n"
	  "implicit none\n"
	  "integer :: ret, f77_test, f90_test\n"
	  "ret = f77_test()\n"
	  "ret = f90_test()\n"
	  "end program test\n")
      endif (lang STREQUAL C)
      # Verify the library works with the current CMAKE compiler
      try_compile(NETCDF_${lang}_FOUND ${scratch_directory} ${test_file}
	CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${NETCDF_INCLUDE_PATH}
	LINK_LIBRARIES ${NETCDF_LIBRARIES}
	OUTPUT_VARIABLE _compile_output)
      unset(scratch_directory)
      unset(test_file)
      unset(_compile_output)
    endif (NOT (NETCDF_INCLUDE_RETURN EQUAL 0) OR
      NOT (NETCDF_LIBRARIES_RETURN EQUAL 0) OR
      NOT (NETCDF_COMPILER_RETURN EQUAL 0) OR
      NOT (NETCDF_COMPILE_FLAGS_RETURN EQUAL 0))

    # Set the return variables
    if (NETCDF_${lang}_FOUND)
      set(NETCDF_${lang}_COMPILER ${NETCDF_COMPILER} CACHE STRING "NetCDF ${lang} compiler" FORCE)
      set(NETCDF_${lang}_COMPILE_FLAGS ${NETCDF_COMPILE_FLAGS} CACHE STRING "NetCDF ${lang} compilation flags" FORCE)
      set(NETCDF_${lang}_INCLUDE_PATH ${NETCDF_INCLUDE_PATH} CACHE STRING "NetCDF ${lang} include path" FORCE)
      set(NETCDF_${lang}_LIBRARIES ${NETCDF_LIBRARIES} CACHE STRING "NetCDF ${lang} linking flags and libraries" FORCE)
      mark_as_advanced(NETCDF_${lang}_COMPILER
	NETCDF_${lang}_COMPILE_FLAGS
	NETCDF_${lang}_INCLUDE_PATH
	NETCDF_${lang}_LIBRARIES)
      set(NetCDF_FIND_REQUIRED_${lang} ${NETCDF_${lang}_FOUND})
    endif (NETCDF_${lang}_FOUND)
  endif (NOT NCDF_CONFIG_${lang})

  # Unset local variables, to clean up namespace
  unset(NETCDF_INCLUDE_RETURN)
  unset(NETCDF_LIBRARIES_RETURN)
  unset(NETCDF_COMPILER_RETURN)
  unset(NETCDF_COMPILE_FLAGS_RETURN)
  unset(NETCDF_COMPILER)
  unset(NETCDF_COMPILE_FLAGS)
  unset(NETCDF_INCLUDE_PATH)
  unset(NETCDF_LIBRARIES)
endforeach (lang C Fortran)

unset(NCDF_CONFIG_C)
unset(NCDF_CONFIG_Fortran)

# Deal with REQUIRED components
foreach (lang C Fortran)
  find_package_handle_standard_args(NetCDF_${lang}
    REQUIRED_VARS NETCDF_${lang}_COMPILER
             NETCDF_${lang}_COMPILE_FLAGS
             NETCDF_${lang}_LIBRARIES
             NETCDF_${lang}_INCLUDE_PATH
    HANDLE_COMPONENTS
    )
endforeach (lang C Fortran)
