#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/mpp directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

skip_test="no"

# Copy file for test.
cp $top_srcdir/test_fms/mpp_io/input_base.nml input.nml


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

if [ $nProc -lt 0 ]
then
   # Couldn't get the number of CPUs, skip the test.
   skip_test="skip"
elif [ $nProc -lt 12 ]
then
   # Need to oversubscribe the MPI
   run_test test_mpp_io 12 $skip_test "true"
fi

run_test test_mpp_io 12 $skip_test
