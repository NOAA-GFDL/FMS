#!/bin/sh

#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# This is part of the GFDL FMS package. This is a shell script to
# execute tests in the test_fms/drifters directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test_common.sh

# Copy file for test_drifters_io test.
cp $top_srcdir/test_fms/drifters/input_base.nml input.nml

# Run tests.

echo "1: Test_drifters_io"
run_test test_drifters_io 2

echo "2: Test_cloud_interpolator"
run_test test_cloud_interpolator 2

echo "3: Test_drifters_comm"
run_test test_drifters_comm 1

echo "4: Test_drifters_core"
run_test test_drifters_core 2

echo "5: Test_quicksort"
run_test test_quicksort 2
