# Precision-based Fortran compiler flags
set(r8_flags "-fdefault-real-8 -fdefault-double-8") # Fortran flags for 64BIT precision
set(r4_flags "-fdefault-real-4") # Fortran flags for 32BIT precision

# GNU Fortran
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcray-pointer -fallow-argument-mismatch -ffree-line-length-none")

set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")

set(CMAKE_Fortran_LINK_FLAGS "" )

#ufs flags
set(CMAKE_Fortran_FLAGS_RELEASEUFS "-O3 -funroll-all-loops -finline-functions")
set(CMAKE_Fortran_FLAGS_DEBUGUFS "-O0 -g -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
