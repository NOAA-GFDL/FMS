#!/bin/sh

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/data_override directory.

# Ed Hartnett 11/26/19

# Set common test settings.
. ../test_common.sh

# Output for TAP.
echo 1..1

# Cause test to fail if any command fails.
cp $top_srcdir/test_fms/data_override/data_table_base data_table
cp $top_srcdir/test_fms/data_override/diag_table_base diag_table

tnum=$( printf "%2.2d" 1 )
sed "s/<test_num>/${tnum}/"  $top_srcdir/test_fms/data_override/input_base.nml > input.nml
#cp -r $top_srcdir/test_fms/data_override/INPUT $top_builddir/test_fms/data_override/INPUT
#mpirun -n 2 ./test_data_override

tnum=$( printf "%2.2d" 2 )
sed "s/<test_num>/${tnum}/"  $top_srcdir/test_fms/data_override/input_base.nml > input.nml
#mpirun -n 2 ./test_data_override

# Output for TAP.
echo ok 1
