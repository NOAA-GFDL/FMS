#!/bin/sh

#***********************************************************************
#*                             Apache License 2.0
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* Licensed under the Apache License, Version 2.0 (the "License");
#* you may not use this file except in compliance with the License.
#* You may obtain a copy of the License at
#*
#*     http://www.apache.org/licenses/LICENSE-2.0
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
#* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#* PARTICULAR PURPOSE. See the License for the specific language
#* governing permissions and limitations under the License.
#***********************************************************************

# Copyright 2021 Seth Underwood

# Set common test settings.
. ../test-lib.sh

# Prepare the directory to run the tests.
touch input.nml

TESTS_SUCCESS='--get-axis-modulo --get-axis-modulo-times --get-axis-cart --lon-in-range --frac-index --nearest-index-increasing --nearest-index-decreasing --axis-edges-increasing --axis_edges-decreasing --tranlon --interp-1d-1d --interp-1d-2d --interp-1d-3d'
TESTS_FAIL='--frac-index-fail --nearest-index-fail'

# TODO: Enable these tests after tranlon's memory corruption bug is fixed.
SKIP_TESTS="test_axis_utils2.19 test_axis_utils2.20"

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
