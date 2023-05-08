#!/bin/sh
# Copyright 2021 Seth Underwood

# Set common test settings.
. ../test-lib.sh

# Prepare the directory to run the tests.
touch input.nml

# Run the test.
test_expect_failure "Test set orbital parameters" '
  mpirun -n 2 ./test_orbital_parameters
'

test_expect_success "Test diurnal_solar" '
  mpirun -n 2 ./test_diurnal_solar
'

test_expect_success "Test daily mean" '
  mpirun -n 2 ./test_daily_mean
'

test_done