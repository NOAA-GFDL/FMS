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
# execute tests in the test_fms/mpp directory.

# Ed Hartnett 11/49/19

# Set common test settings.
. ../test-lib.sh

# Copy file for test.
touch input.nml

rm -rf INPUT
mkdir INPUT


# Run the test.
test_expect_success "test time interpolation with r8_kind" '
  mpirun -n 4 ./test_time_interp_r8
'
test_expect_success "test time interpolation with r4_kind" '
  mpirun -n 4 ./test_time_interp_r4
'

rm -rf INPUT
mkdir INPUT

# nml for calender type
cat <<_EOF > input.nml
&test_time_interp_external_nml
cal_type="julian"
/
_EOF

test_expect_success "test time interpolation external with r8_kind (julian)" '
  mpirun -n 4 ./test_time_interp_external_r8
'

test_expect_success "test time interpolation external with r4_kind (julian)" '
  mpirun -n 4 ./test_time_interp_external_r4
'
sed -i 's/julian/no_leap/' input.nml

test_expect_success "test time interpolation external with r8_kind (no_leap)" '
  mpirun -n 4 ./test_time_interp_external_r8
'
test_expect_success "test time interpolation external with r4_kind (no_leap)" '
  mpirun -n 4 ./test_time_interp_external_r4
'

test_expect_success "test time interpolation external with conservative horizontal interp r4_kind" '
  mpirun -n 4 ./test_time_interp_conservative_hi_r4
'

test_expect_success "test time interpolation external with conservative horizontal interp r8_kind" '
  mpirun -n 4 ./test_time_interp_conservative_hi_r8
'


test_done
