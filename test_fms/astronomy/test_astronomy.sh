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
# execute tests in the test_fms/astronomy directory.

# Caitlyn McAllister

# Set common test settings.
. ../test-lib.sh

# Prepare the directory to run the tests.
touch input.nml

# Run the test.

test_expect_failure "Test set orbital parameters - test 1: r4_kind" '
  mpirun -n 1 ./test_set_orbital_parameters_r4
'

test_expect_failure "Test set orbital parameters - test 1: r8_kind" '
  mpirun -n 1 ./test_set_orbital_parameters_r8
'

test_expect_success "Test get orbital parameters - test 2: r4_kind" '
  mpirun -n 1 ./test_get_orbital_parameters_r4
'

test_expect_success "Test get orbital parameters - test 2: r8_kind" '
  mpirun -n 1 ./test_get_orbital_parameters_r8
'

test_expect_success "Test diurnal_solar 2d-0d - test 3: r4_kind" '
  mpirun -n 1 ./test_diurnal_solar_r4
'

test_expect_success "Test diurnal_solar 2d-0d - test 3: r8_kind" '
  mpirun -n 1 ./test_diurnal_solar_r8
'

test_expect_success "Test optional argument dt for diurnal_solar 2d-0d - test 4: r4_kind" '
  mpirun -n 1 ./test_dt_option_diurnal_r4
'

test_expect_success "Test optional argument dt for diurnal_solar 2d-0d - test 4: r8_kind" '
  mpirun -n 1 ./test_dt_option_diurnal_r8
'

test_expect_success "Test diurnal_solar_cal 2d-0d - test 5: r4_kind" '
  mpirun -n 1 ./test_diurnal_solar_cal_r4
'

test_expect_success "Test diurnal_solar_cal 2d-0d - test 5: r8_kind" '
  mpirun -n 1 ./test_diurnal_solar_cal_r8
'

test_expect_success "Test daily mean 2d-0d - test 6: r4_kind" '
  mpirun -n 1 ./test_daily_mean_r4
'

test_expect_success "Test daily mean 2d-0d - test 6: r8_kind" '
  mpirun -n 1 ./test_daily_mean_r8
'

test_expect_success "Test daily mean cal 2d-0d - test 7: r4_kind" '
  mpirun -n 1 ./test_daily_mean_cal_r4
'

test_expect_success "Test daily mean cal 2d-0d - test 7: r8_kind" '
  mpirun -n 1 ./test_daily_mean_cal_r8
'

test_expect_success "Test annual mean 2d - test 8: r4_kind" '
  mpirun -n 1 ./test_annual_mean_solar_2d_r4
'

test_expect_success "Test annual mean 2d - test 8: r8_kind" '
  mpirun -n 1 ./test_annual_mean_solar_2d_r8
'

test_expect_success "Test annual mean 1d - test 9: r4_kind" '
  mpirun -n 2 ./test_annual_mean_solar_1d_r4
'

test_expect_success "Test annual mean 1d - test 9: r8_kind" '
  mpirun -n 2 ./test_annual_mean_solar_1d_r8
'

test_expect_success "Test annual mean 2level - test 10: r4_kind" '
  mpirun -n 2 ./test_annual_mean_solar_2level_r4
'

test_expect_success "Test annual mean 2level - test 10: r8_kind" '
  mpirun -n 2 ./test_annual_mean_solar_2level_r8
'

test_done