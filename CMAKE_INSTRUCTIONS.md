# Instructions for building FMS with cmake
##############################
# 1 Environment Variables
# Set the environment variables FC, CC, FFLAGS, and CPPFLAGS
##############################
# Use FMS Autotools files for generating configuration information for cmake:

$ cd FMS
$ autoreconf -i
$ ./configure

Notes: The above commands create the files FMSConfig.cmake.in and FMSConfigVersion.cmake.in. Also, the autotools command "make" and "make install" are not run.
##########################################
# Have cmake build and install FMS.
# (Below, fms_install_path is the full install directory path name
# you select):

$ cd FMS
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX:PATH=fms_install_path ..
$ cmake --build . --target install --config Release

# When the above command finishes, the fms_install_path will have an include
# and a lib directory. The lib directory will have these files:
#     libFMS.a
#     /cmake/FMS/FMSTargets.cmake
#     /cmake/FMS/FMSTargets-noconfig.cmake
#     /cmake/FMS/FMSConfig.cmake
#     /cmake/FMS/FMSConfigVersion.cmake
######################
