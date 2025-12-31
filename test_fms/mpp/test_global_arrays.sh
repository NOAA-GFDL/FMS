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

# Ryan Mulhall 2/2021

# Set common test settings.
. ../test-lib.sh

cat <<_EOF  > input.nml
&test_global_arrays_nml
  test_simple = .true.
  test_full = .false.
/
_EOF

test_expect_success "mpp_global_sum/max/min with simple domain" '
    mpirun -n 8 ./test_global_arrays
'

cat <<_EOF  > input.nml
&test_global_arrays_nml
  test_simple = .false.
  test_full = .true.
/
_EOF

test_expect_success "mpp_global_sum/max/min with symmetry and halos" '
    mpirun -n 6 ./test_global_arrays
'

test_done
