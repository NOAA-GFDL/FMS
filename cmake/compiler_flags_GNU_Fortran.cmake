# Precision-based Fortran compiler flags
set(r8_flags "-fdefault-real-8 -fdefault-double-8") # Fortran flags for 64BIT precision

# GNU Fortan
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcray-pointer -fconvert=big-endian -ffree-line-length-none -fno-range-check -fbacktrace")

set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops -finline-functions")

set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )

set(CMAKE_Fortran_LINK_FLAGS "" )
