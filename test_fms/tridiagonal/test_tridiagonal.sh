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
# execute tests in the test_fms/time_manager directory.

# Ryan Mulhall 9/2023

# Set common test settings.
. ../test-lib.sh

rm -f input.nml && touch input.nml

test_expect_success "test tridiagonal functionality 32 bit reals" '
    mpirun -n 1 ./test_tridiagonal_r4
'
test_expect_success "test tridiagonal functionality 64 bit reals" '
    mpirun -n 1 ./test_tridiagonal_r8
'
# tries to call without a,b,c args provided or previously set
cat <<_EOF > input.nml
&test_tridiagonal_nml
do_error_check = .true.
/
_EOF
test_expect_failure "error out if passed in incorrect real size (r4_kind)" '
    mpirun -n 1 ./test_tridiagonal_r4
'
test_expect_failure "error out if passed in incorrect real size (r8_kind)" '
    mpirun -n 1 ./test_tridiagonal_r8
'

test_done
