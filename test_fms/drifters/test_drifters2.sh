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
