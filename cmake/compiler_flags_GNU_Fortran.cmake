# Precision-based Fortran compiler flags
set(r8_flags "-fdefault-real-8 -fdefault-double-8") # Fortran flags for 64BIT precision

# GNU Fortran
set(CMAKE_Fortran_FLAGS "")

set(base_flags "-fcray-pointer -fdefault-double-8 -Waliasing -ffree-line-length-none -fno-range-check -fallow-argument-mismatch")

set(CMAKE_Fortran_FLAGS_RELEASE "${base_flags} -O2 -fno-expensive-optimizations")
set(CMAKE_Fortran_FLAGS_DEBUG "${base_flags} -O0 -g -W -fbounds-check -ffpe-trap=invalid,zero,overflow")

set(CMAKE_Fortran_LINK_FLAGS "" )

#ufs flags
set(CMAKE_Fortran_FLAGS_RELEASEUFS "-O3 -funroll-all-loops -finline-functions")
set(CMAKE_Fortran_FLAGS_DEBUGUFS "-O0 -g -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
