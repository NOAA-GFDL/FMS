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
# execute tests in the test_fms/time_manager directory.

# Ed Hartnett 11/29/19

# Set common test settings.
. ../test-lib.sh

rm input.nml && touch input.nml

test_expect_success "test tridiagonal functionality 32 bit reals" '
    mpirun -n 1 ./test_tridiagonal_r4
'
test_expect_success "test tridiagonal functionality 64 bit reals" '
    mpirun -n 1 ./test_tridiagonal_r8
'
# tries to call without a,b,c args provided or preciously set
cat <<_EOF > input.nml
&test_tridiagonal_nml
do_error_check = .true.
/
_EOF
test_expect_failure "test error out on incorrect kind input" '
    mpirun -n 1 ./test_tridiagonal_r8
'
test_expect_failure "test error out on incorrect kind input" '
    mpirun -n 1 ./test_tridiagonal_r4
'

test_done
