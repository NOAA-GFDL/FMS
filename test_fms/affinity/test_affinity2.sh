#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/axis_utils directory.

# Ed Hartnett 11/26/19

# Set common test settings.
. ../test_common.sh

# Copy and rename namelist file.
cp $top_srcdir/test_fms/affinity/input_base.nml input.nml

# Only run the test if on slurm:
if test "x$mpi_launcher" == "xsrun" ; then
   run_test test_affinity 6
fi
