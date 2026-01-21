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
# execute tests in the test_fms/drifters directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test-lib.sh

# input.nml for tests
cat <<_EOF > input.nml
&drifters_comm_nml
   nx = 11
  ny = 21
  halox = 2
  haloy = 2
  u0    = 1.0
  v0    = 0.0
  dt    = 0.1
  nt    = 10
/
_EOF

# Run tests.
test_expect_success "Test drifters IO" '
  mpirun -n 2 ./test_drifters_io
'

test_expect_success "Test cloud interpolator" '
  mpirun -n 2 ./test_cloud_interpolator
'

test_expect_success "Test drifters comm" '
  mpirun -n 1 ./test_drifters_comm
'

test_expect_success "Test drifters core" '
  mpirun -n 2 ./test_drifters_core
'

test_expect_success "test quicksort" '
  mpirun -n 2 ./test_quicksort
'

test_done
