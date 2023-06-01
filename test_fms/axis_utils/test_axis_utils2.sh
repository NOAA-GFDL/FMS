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

# Copyright 2021 Seth Underwood

# Set common test settings.
. ../test-lib.sh

# Prepare the directory to run the tests.
touch input.nml

TESTS_SUCCESS='--get-axis-modulo --get-axis-modulo-times --get-axis-cart --lon-in-range --frac-index --nearest-index-increasing --nearest-index-decreasing --axis-edges --tranlon --interp-1d-1d --interp-1d-2d --interp-1d-3d'
TESTS_FAIL='--frac-index-fail --nearest-index-fail'

# TODO: Enable these tests after tranlon's memory corruption bug is fixed.
SKIP_TESTS="test_axis_utils2.15 test_axis_utils2.16"

# Run the tests

for t in $TESTS_SUCCESS
do
  r4cmd="./test_axis_utils_r4 $t"
  r8cmd="./test_axis_utils_r8 $t"

  test_expect_success "Testing axis utils: $r4cmd" "mpirun -n 1 $r4cmd"
  test_expect_success "Testing axis utils: $r8cmd" "mpirun -n 1 $r8cmd"
done

for t in $TESTS_FAIL
do
  r4cmd="./test_axis_utils_r4 $t"
  r8cmd="./test_axis_utils_r8 $t"

  test_expect_failure "Testing axis utils: $r4cmd" "mpirun -n 1 $r4cmd"
  test_expect_failure "Testing axis utils: $r8cmd" "mpirun -n 1 $r8cmd"
done

test_done
