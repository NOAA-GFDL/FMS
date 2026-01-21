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

test_expect_success "test time interpolation external with conservative horizontal interp r8_kind" '
  mpirun -n 4 ./test_time_interp_conservative_hi_r8
'
test_expect_success "test time interpolation external with conservative horizontal interp r4_kind" '
  mpirun -n 4 ./test_time_interp_conservative_hi_r4
'

test_done
