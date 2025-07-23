# Precision-based Fortran compiler flags
set(r8_flags "-r8") # Fortran flags for 64BIT precision

# GNU Fortan
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ")

set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast")

set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )

set(CMAKE_Fortran_LINK_FLAGS "" )
