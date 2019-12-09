#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/drifters directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# Copy file for test_drifters_io test.
cp $top_srcdir/test_fms/drifters/input_base.nml input.nml

# Run tests.

if test "$mpi_launcher" != "" ; then 
    npes="-n 2"
fi

echo "1: Test_drifters_io"
$mpi_launcher $npes ./test_drifters_io

echo "2: Test_cloud_interpolator"
$mpi_launcher $npes ./test_cloud_interpolator

echo "3: Test_drifters_comm"
./test_drifters_comm

echo "4: Test_drifters_core"
$mpi_launcher $npes ./test_drifters_core

echo "5: Test_quicksort"
echo $mpi_launcher $npes ./test_quicksort
