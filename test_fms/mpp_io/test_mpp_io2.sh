#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/mpp directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# Copy file for test.
cp $top_srcdir/test_fms/mpp_io/input_base.nml input.nml

# If there is a mpi launcher set the number of processors otherwise just ./
if test "$mpi_launcher" != "" ; then 
# Get the number of available CPUs on the system
    if [ $(command -v nproc) ]
    then
    # Looks like a linux system
        nProc=$(nproc)
    elif [ $(command -v sysctl) ]
    then
        # Looks like a Mac OS X system
        nProc=$(sysctl -n hw.physicalcpu)
    else
        nProc=-1
    fi
    #Try to set the numbers of processors to 12 if possible (skylake needs 12 to pass ...)
    if [ $nProc -gt 12 ]
    then  
        npes="-n 12"
    fi
fi

$mpi_launcher $npes ./test_mpp_io
