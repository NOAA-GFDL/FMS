# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# compile C with /opt/oneapi/2022.2/mpi/2021.6.0/bin/mpiicc
# compile Fortran with /opt/oneapi/2022.2/mpi/2021.6.0/bin/mpiifort
C_DEFINES = -DENABLE_QUAD_PRECISION -DGFDL_CONSTANTS -DINTERNAL_FILE_NML -DOVERLOAD_R4 -DOVERLOAD_R8 -Duse_libMPI -Duse_netCDF

C_INCLUDES = -I/home/Caitlyn.Mcallister/axisUtils/include -I/home/Caitlyn.Mcallister/axisUtils/fms -I/home/Caitlyn.Mcallister/axisUtils/fms2_io/include -I/home/Caitlyn.Mcallister/axisUtils/mpp/include -I/home/Caitlyn.Mcallister/axisUtils/axis_utils

C_FLAGS = -I/opt/netcdf/4.8.0/ONEAPI/include -I/opt/hdf5/1.12.0/ONEAPI/include -sox -traceback -qno-opt-dynamic-align -O2 -debug minimal -fp-model source

Fortran_DEFINES = -DENABLE_QUAD_PRECISION -DGFDL_CONSTANTS -DINTERNAL_FILE_NML -DOVERLOAD_R4 -DOVERLOAD_R8 -Duse_libMPI -Duse_netCDF

Fortran_INCLUDES = -I/home/Caitlyn.Mcallister/axisUtils/include -I/home/Caitlyn.Mcallister/axisUtils/fms -I/home/Caitlyn.Mcallister/axisUtils/fms2_io/include -I/home/Caitlyn.Mcallister/axisUtils/mpp/include -I/home/Caitlyn.Mcallister/axisUtils/axis_utils -I/opt/netcdf/4.8.0/ONEAPI/include

Fortran_FLAGS =  -fpp -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -align array64byte -nowarn -sox -traceback -O2 -debug minimal -fp-model source -nowarn -qoverride-limits -qno-opt-dynamic-align -qopt-prefetch=3

