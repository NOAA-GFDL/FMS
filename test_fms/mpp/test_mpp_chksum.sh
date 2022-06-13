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

echo "&test_mpp_chksum_nml" > input.nml
echo "test_num = 1" >> input.nml
# replaces defaults with smaller sizes if stack size is limited
if [ $STACK_LIMITED ]; then
  echo "nx = 64" >> input.nml
  echo "ny = 64" >> input.nml
fi
echo "/" >> input.nml

test_expect_success "mpp_chksum simple functionality" '
    mpirun -n 4 ./test_mpp_chksum
'

sed -i 's/test_num = 1/test_num = 2/' input.nml
test_expect_success "mpp integer checksums" '
    mpirun -n 4 ./test_mpp_chksum
'

sed -i 's/test_num = 2/test_num = 3/' input.nml
test_expect_success "mpp_chksum with mixed precision" '
    mpirun -n 4 ./test_mpp_chksum
'
test_done
