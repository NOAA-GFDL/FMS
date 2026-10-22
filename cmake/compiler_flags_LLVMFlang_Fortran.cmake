# Precision-based Fortran compiler flags
set(r8_flags "-fdefault-real-8 -fdefault-double-8")
set(r4_flags "-fdefault-real-4")

# LLVM Flang
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
set(CMAKE_Fortran_LINK_FLAGS "")

# UFS flags are not known for LLVM Flang.
