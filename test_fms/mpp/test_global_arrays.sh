#!/bin/sh

#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the GFDL Flexible Modeling System (FMS).
#
# FMS is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FMS is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
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
