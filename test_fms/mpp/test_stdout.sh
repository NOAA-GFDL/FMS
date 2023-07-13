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

# Colin Gladue 06/12/20
# Ryan Mulhall 02/2021

# Set common test settings.
. ../test-lib.sh

# ensure input.nml file present
touch input.nml

# Run test with one processor
test_expect_success "get stdout with 1 PE" '
    mpirun -n 2 ./test_stdout
'
test_expect_success "get stdout with 2 PEs" '
    mpirun -n 2 ./test_stdout
'
test_done
