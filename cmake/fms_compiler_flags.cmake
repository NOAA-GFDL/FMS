################################################################################
# Fortran
################################################################################

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  include(compiler_flags_GNU_Fortran)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  include(compiler_flags_Intel_Fortran)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  include(compiler_flags_IntelLLVM_Fortran)
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

################################################################################
# C
################################################################################

if(CMAKE_C_COMPILER_ID MATCHES "GNU")
  include(compiler_flags_GNU_C)
elseif(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  include(compiler_flags_Intel_C)
elseif(CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
  include(compiler_flags_IntelLLVM_C)
elseif(CMAKE_C_COMPILER_ID MATCHES "Clang")
  include(compiler_flags_Clang_C)
else()
  message(WARNING "C compiler with ID ${CMAKE_C_COMPILER_ID} will be used with CMake default options")
endif()
